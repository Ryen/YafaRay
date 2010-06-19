#include <integrators/SPPM.h>
#include <sstream>
#include <cmath>
#include <algorithm>
__BEGIN_YAFRAY


SPPM::SPPM(unsigned int dPhotons, int _passnum)
{
	type = SURFACE;
	intpb = 0;
	integratorName = "SPPM";
	integratorShortName = "SPPM";
	nPhotons = dPhotons;
	//nSearch = dPhotons; // need to find all photons in the region. But this setting may cause stack overflow;
	passNum = _passnum;
	totalnPhotons = 0;
	initialFactor = 1.f;
}

SPPM::~SPPM()
{}

bool SPPM::preprocess()
{
	return true;
}
bool SPPM::render(yafaray::imageFilm_t *image)
{
	std::stringstream passString;
	imageFilm = image;

	//scene->getAAParameters(AA_samples, AA_passes, AA_inc_samples, AA_threshold);
	//iAA_passes = 1.f / (float) AA_passes;
	//Y_INFO << integratorName << ": Rendering " << AA_passes << " passes" << std::endl;
	//Y_INFO << integratorName << ": Min. " << AA_samples << " samples" << std::endl;
	//Y_INFO << integratorName << ": "<< AA_inc_samples << " per additional pass" << std::endl;
	//Y_INFO << integratorName << ": Max. " << AA_samples + std::max(0,AA_passes-1) * AA_inc_samples << " total samples" << std::endl;
	//passString << "Rendering pass 1 of " << std::max(1, AA_passes) << "...";
	//Y_INFO << integratorName << ": " << passString.str() << std::endl;

	//passString << "Render Start!!!" <<std::endl;

	//if(intpb) intpb->setTag("Render Start!!!"); 

	gTimer.addEvent("rendert");
	gTimer.start("rendert");
	imageFilm->init(passNum); // need to doule-check how it effect sppm.
	
	const camera_t* camera = scene->getCamera();
	progressiveData.reserve(camera->resX() * camera->resY()); // used for SPPM

	maxDepth = 0.f;
	minDepth = 1e38f;

	point3d_t a(1e38f,1e38f,1e38f);
	point3d_t g(0.f, 0.f, 0.f); // Is it right for the initial g point?
	bound_t bBox(a, g);

	if(scene->doDepth()) //  Use it to set get scene BBox.
	{
		const camera_t* camera = scene->getCamera();
		diffRay_t c_ray;
		int end_x=camera->resX(), end_y=camera->resY();
		float wt = 0.f;
		surfacePoint_t sp;

		for(int i=0; i<end_y; ++i)
		{
			for(int j=0; j<end_x; ++j)
			{
				c_ray.tmax = -1.f;
				c_ray = camera->shootRay(i, j, 0.5f, 0.5f, wt);
				scene->intersect(c_ray, sp);
				if(c_ray.tmax > maxDepth) maxDepth = c_ray.tmax;
				if(c_ray.tmax < minDepth && c_ray.tmax >= 0.f) minDepth = c_ray.tmax;
				bBox.include(sp.P);
			}
		}
		if(maxDepth > 0.f) maxDepth = 1.f / (maxDepth - minDepth);
	}

	bBox = scene->getSceneBound(); // Now using Scene Bound 
	// initialize SPPM statistics
	float initialRadius = ((bBox.longX() + bBox.longY() + bBox.longZ()) / 3.f) / ((camera->resX() + camera->resY()) / 2.0f) * 2.f ;

	initialRadius = std::min(initialRadius, 1.f); //Fix the overflow bug

	int size = camera->resX() * camera->resY();
	for(int i = 0; i < size; i++)
	{
		shared_statistics ts;

		ts.accflux = colorA_t(0.f);
		ts.photoncount = 0;
		ts.radius2 = (initialRadius * initialFactor) * (initialRadius * initialFactor);

		ts.acclightsource = colorA_t(0.f);
		ts.constanthits = 0;
		ts.surfacehits = 0;

		progressiveData.push_back(ts);
	}

	renderPass(1, 0, false);

	for(int i=0; i<passNum; ++i) // progress pass
	{
		if(scene->getSignals() & Y_SIG_ABORT) break;
		//imageFilm->setAAThreshold(AA_threshold);
		imageFilm->nextPass(false, integratorName);
		nRefined = 0;
		renderPass(1, 1 + (i-1)*1, false); // offset seems only used for sampling?
		std::cout<<  "This pass refined "<<nRefined<<" of "<<size<<" pixels."<<"\n";
	}
	maxDepth = 0.f;
	gTimer.stop("rendert");
	Y_INFO << integratorName << ": Overall rendertime: "<< gTimer.getTime("rendert")<<"s\n";

	return true;
}

bool SPPM::renderPass(int samples,int offset, bool adaptive)
{
	prePass(samples, offset, adaptive);
	renderArea_t a;
	while(imageFilm->nextArea(a))
	{
		if(scene->getSignals() & Y_SIG_ABORT) break;
		renderTile(a, samples, offset, adaptive,0);
		imageFilm->finishArea(a);
	}
	return true;
}

bool SPPM::renderTile(renderArea_t &a, int n_samples, int offset, bool adaptive, int threadID)
{
	int x, y;
	const camera_t* camera = scene->getCamera();
	bool do_depth = scene->doDepth();
	x=camera->resX();
	y=camera->resY();
	diffRay_t c_ray;
	ray_t d_ray;
	PFLOAT dx=0.5, dy=0.5, d1=1.0/(PFLOAT)n_samples;
	float lens_u=0.5f, lens_v=0.5f;
	PFLOAT wt, wt_dummy;
	random_t prng(offset*(x*a.Y+a.X)+123);
	renderState_t rstate(&prng);
	rstate.threadID = threadID;
	rstate.cam = camera;
	bool sampleLns = camera->sampleLense();
	int pass_offs=offset, end_x=a.X+a.W, end_y=a.Y+a.H;
	
	for(int i=a.Y; i<end_y; ++i)
	{
		for(int j=a.X; j<end_x; ++j)
		{
			if(scene->getSignals() & Y_SIG_ABORT) break;

			rstate.pixelNumber = x*i+j;
			rstate.samplingOffs = fnv_32a_buf(i*fnv_32a_buf(j));//fnv_32a_buf(rstate.pixelNumber);
			float toff = scrHalton(5, pass_offs+rstate.samplingOffs); // **shall be just the pass number...**

			for(int sample=0; sample<n_samples; ++sample) //set n_samples = 1.
			{
				rstate.setDefaults();
				rstate.pixelSample = pass_offs+sample;
				rstate.time = addMod1((PFLOAT)sample*d1, toff); //(0.5+(PFLOAT)sample)*d1;
				// the (1/n, Larcher&Pillichshammer-Seq.) only gives good coverage when total sample count is known
				// hence we use scrambled (Sobol, van-der-Corput) for multipass AA

				dx = RI_S(rstate.pixelSample, rstate.samplingOffs);
				dy = RI_vdC(rstate.pixelSample, rstate.samplingOffs);

				if(sampleLns)
				{
					lens_u = scrHalton(3, rstate.pixelSample+rstate.samplingOffs);
					lens_v = scrHalton(4, rstate.pixelSample+rstate.samplingOffs);
				}
				c_ray = camera->shootRay(j+dx, i+dy, lens_u, lens_v, wt); // wt need to be considered
				if(wt==0.0)
				{
					imageFilm->addSample(colorA_t(0.f), j, i, dx, dy, &a); //maybe not need
					continue;
				}
				//setup ray differentials
				d_ray = camera->shootRay(j+1+dx, i+dy, lens_u, lens_v, wt_dummy);
				c_ray.xfrom = d_ray.from;
				c_ray.xdir = d_ray.dir;
				d_ray = camera->shootRay(j+dx, i+1+dy, lens_u, lens_v, wt_dummy);
				c_ray.yfrom = d_ray.from;
				c_ray.ydir = d_ray.dir;
				c_ray.time = rstate.time;
				c_ray.hasDifferentials = true;
				// col = T * L_o + L_v
				diffRay_t c_ray_copy = c_ray;

				//for sppm progressive
				int index = i*camera->resX() + j; 
				shared_statistics &temp = progressiveData[index];

				curPhotons = 0;
				constantColor = colorA_t(0.f);
				curRadius2 = temp.radius2;
				colorA_t flux = integrate(rstate, c_ray); // L_o
				flux *= scene->volIntegrator->transmittance(rstate, c_ray); // T
				flux += scene->volIntegrator->integrate(rstate, c_ray); // L_v

				// progressive refinement
				const float _alpha = 0.7;
				float g = 1.0f;

				// The author's refine formular
				//if(curPhotons != 0)
				//{
				//	float g = std::min((temp.photoncount + _alpha * curPhotons) / (temp.photoncount + curPhotons), 1.0f);
				//	temp.radius = temp.radius * sqrt(g);
				//	temp.photoncount += curPhotons * _alpha;
				//	temp.accflux = (temp.accflux + flux) * g;		
				//}

				//Dade's formula, need to fix
				if(curPhotons > 0)
				{
						const unsigned long long pcount = temp.photoncount + curPhotons;
						const float g = _alpha * pcount / (temp.photoncount * _alpha + curPhotons);
						temp.photoncount = pcount;
						temp.accflux = (temp.accflux + flux) * g;
						temp.radius2 *= g;
						nRefined++;
				}
				//radiance estimate
				colorA_t color = temp.accflux / (temp.radius2 * M_PI * totalnPhotons);
				imageFilm->addSample(wt * (color+constantColor), j, i, dx, dy, &a);

				if(do_depth) // need to know how it works
				{
					float depth = 0.f;
					if(c_ray.tmax > 0.f)
					{
						depth = 1.f - (c_ray.tmax - minDepth) * maxDepth; // Distance normalization
					}
					imageFilm->addDepthSample(0, depth, j, i, dx, dy);
				}
			}
		}
	}
	return true;
}

//photon pass, scatter photon 
void SPPM::prePass(int samples, int offset, bool adaptive)
{
	std::stringstream set;
	gTimer.addEvent("prePass");
	gTimer.start("prePass");

	Y_INFO << integratorName << ": Starting Photon tracing pass...\n";

	//if(trShad)
	//{
	//	set << "ShadowDepth [" << sDepth << "]";
	//}
	//if(!set.str().empty()) set << "+";
	//set << "RayDepth [" << rDepth << "]";

	diffuseMap.clear();
	background = scene->getBackground();
	lights = scene->lights;
	std::vector<light_t*> tmplights;

	//if(background)
	//{
	//	light_t *bgl = background->getLight();
	//	if(bgl)
	//	{
	//		lights.push_back(bgl);
	//		hasBGLight = true;
	//		if(!set.str().empty()) set << "+";
	//		set << "IBL";
	//	}
	//}
	
	settings = set.str();
	
	ray_t ray;
	float lightNumPdf, lightPdf, s1, s2, s3, s4, s5, s6, s7, sL;
	int numCLights = 0;
	int numDLights = 0;
	float fNumLights = 0.f;
	float *energies = NULL;
	color_t pcol;

	tmplights.clear();

	for(int i=0;i<(int)lights.size();++i)
	{
		numDLights++;
		tmplights.push_back(lights[i]);
	}
	
	fNumLights = (float)numDLights;
	energies = new float[numDLights];

	for(int i=0;i<numDLights;++i) energies[i] = tmplights[i]->totalEnergy().energy();

	lightPowerD = new pdf1D_t(energies, numDLights);
	
	Y_INFO << integratorName << ": Light(s) photon color testing for photon map:\n";
	for(int i=0;i<numDLights;++i)
	{
		pcol = tmplights[i]->emitPhoton(.5, .5, .5, .5, ray, lightPdf);
		lightNumPdf = lightPowerD->func[i] * lightPowerD->invIntegral;
		pcol *= fNumLights*lightPdf/lightNumPdf; //remember that lightPdf is the inverse of the pdf, hence *=...
		Y_INFO << integratorName << ": Light ["<<i+1<<"] Photon col:"<<pcol<<" | lnpdf: "<<lightNumPdf<<"\n";
	}
	
	delete[] energies;
	
	//shoot photons
	bool done=false;
	unsigned int curr=0;

	
	surfacePoint_t sp;
	renderState_t state;
	unsigned char userdata[USER_DATA_SIZE+7];
	state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
	state.cam = scene->getCamera();
	progressBar_t *pb;
	int pbStep;
	if(intpb) 
		pb = intpb;
	else 
		pb = new ConsoleProgressBar_t(80);
	
	Y_INFO << integratorName << ": Building photon map...\n";
	
	pb->init(128);
	pbStep = std::max(1U, nPhotons/128);
	pb->setTag("Building photon map...");
	//Pregather diffuse photons

	float invDiffPhotons = 1.f / (float)nPhotons;
	
	while(!done)
	{
		if(scene->getSignals() & Y_SIG_ABORT) {  pb->done(); if(!intpb) delete pb; return; }
		//state.chromatic = true;
		//state.wavelength = scrHalton(5, curr);
		
		//s1 = RI_vdC(curr+offset);
		//s2 = scrHalton(2, curr+offset);
		//s3 = scrHalton(3, curr+offset);
		//s4 = scrHalton(4, curr+offset);

		s1 = ourRandom();
		s2 = ourRandom();
		s3 = ourRandom();
		s4 = ourRandom();


		sL = float(curr) * invDiffPhotons;
		int lightNum = lightPowerD->DSample(sL, &lightNumPdf);
		if(lightNum >= numDLights){ Y_ERROR << integratorName << ": lightPDF sample error! "<<sL<<"/"<<lightNum<<"... stopping now.\n"; delete lightPowerD; return; }

		pcol = tmplights[lightNum]->emitPhoton(s1, s2, s3, s4, ray, lightPdf);
		ray.tmin = MIN_RAYDIST;
		ray.tmax = -1.0;
		pcol *= fNumLights*lightPdf/lightNumPdf; //remember that lightPdf is the inverse of th pdf, hence *=...

		if(pcol.isBlack())
		{
			++curr;
			done = (curr >= nPhotons);
			continue;
		}

		int nBounces=0;
		bool causticPhoton = false;
		bool directPhoton = true;
		const material_t *material = NULL;
		BSDF_t bsdfs;

		while( scene->intersect(ray, sp) ) //scatter photons.
		{
			if(isnan(pcol.R) || isnan(pcol.G) || isnan(pcol.B))
			{ Y_WARNING << integratorName << ": NaN  on photon color for light" << lightNum + 1 << ".\n"; continue; }
			
			color_t transm(1.f);
			color_t vcol(0.f);
			const volumeHandler_t* vol;
			
			if(material)
			{
				if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * -ray.dir < 0)))
				{
					if(vol->transmittance(state, ray, vcol)) transm = vcol;
				}
			}
			
			vector3d_t wi = -ray.dir, wo;
			material = sp.material;
			material->initBSDF(state, sp, bsdfs);
			
			if(bsdfs & (BSDF_DIFFUSE))
			{
				//deposit photon on surface
				if( 1 || !causticPhoton)
				{
					photon_t np(wi, sp.P, pcol);// pcol used here
					np.shadeN = sp.Ng;
					diffuseMap.pushPhoton(np);
					diffuseMap.setNumPaths(curr); // fixed it
				}

			}
			// need to break in the middle otherwise we scatter the photon and then discard it => redundant
			//if(nBounces == maxBounces)break;  

			// scatter photon
			//int d5 = 3*nBounces + 5;

			//s5 = scrHalton(d5, curr);
			//s6 = scrHalton(d5+1, curr);
			//s7 = scrHalton(d5+2, curr);

			s5 = ourRandom();
			s6 = ourRandom();
			s7 = ourRandom();
			
			pSample_t sample(s5, s6, s7, BSDF_ALL, pcol, transm);

			bool scattered = material->scatterPhoton(state, sp, wi, wo, sample);
			if(!scattered) break; //photon was absorped.  

			pcol = sample.color;

			causticPhoton = ((sample.sampledFlags & (BSDF_GLOSSY | BSDF_SPECULAR | BSDF_DISPERSIVE)) && directPhoton) ||
							((sample.sampledFlags & (BSDF_GLOSSY | BSDF_SPECULAR | BSDF_FILTER | BSDF_DISPERSIVE)) && causticPhoton);
			directPhoton = (sample.sampledFlags & BSDF_FILTER) && directPhoton;

			ray.from = sp.P;
			ray.dir = wo;
			ray.tmin = MIN_RAYDIST;
			ray.tmax = -1.0;
			++nBounces;

			// Russian Roulette
			if(nBounces > 3)
			{
				float continueProbability = 0.75f;
				if (ourRandom() > continueProbability) break;
				pcol *= 1.f / continueProbability;
			}

		}
		++curr;
		if(curr % pbStep == 0) pb->update();
		done = (curr >= nPhotons);
	}
	pb->done();
	pb->setTag("Photon map built.");
	Y_INFO << integratorName << ":Photon map built.\n";
	Y_INFO << integratorName << ": Shot "<<curr<<" photons from " << numDLights << " light(s)\n";

	delete lightPowerD;

	totalnPhotons +=  nPhotons; // used for SPPM

	delete[] energies;

	Y_INFO << integratorName << ": Stored photons: "<<diffuseMap.nPhotons()<<"\n";
	
	if(diffuseMap.nPhotons() > 0) //
	{
		Y_INFO << integratorName << ": Building photons kd-tree:\n";
		pb->setTag("Building photons kd-tree...");
		diffuseMap.updateTree();
		Y_INFO << integratorName << ": Done.\n";
	}

	if(diffuseMap.nPhotons() < 50) { Y_ERROR << integratorName << ": Too few photons, stopping now.\n"; return; }
	tmplights.clear();

	if(!intpb) delete pb;
	
	gTimer.stop("prePass");
	Y_INFO << integratorName << ": Photonmap building time: " << gTimer.getTime("prePass") << "\n";

	return;
}

//collect photon information, it should return non-normalized flux for each pixel.
colorA_t SPPM::integrate(renderState_t &state, diffRay_t &ray/*, sampler_t &sam*/) const
{
	static int _nMax=0;
	static int calls=0;
	++calls;
	color_t col(0.0);
	CFLOAT alpha=0.0;
	surfacePoint_t sp;
	
	void *o_udat = state.userdata;
	bool oldIncludeLights = state.includeLights;
	if(scene->intersect(ray, sp))
	{
		unsigned char userdata[USER_DATA_SIZE+7];
		state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
		if(state.raylevel == 0)
		{
			state.chromatic = true;
			state.includeLights = true;
		}
		BSDF_t bsdfs;
		vector3d_t N_nobump = sp.N;
		vector3d_t wo = -ray.dir;
		const material_t *material = sp.material;
		material->initBSDF(state, sp, bsdfs);
		constantColor += material->emit(state, sp, wo);
		state.includeLights = false;
		spDifferentials_t spDiff(sp, ray);

        // remove FG here
		//if(bsdfs & BSDF_DIFFUSE) directColor += estimateDirect_PH(state, sp, lights, scene, wo, false, sDepth);
		
		if(0)
		{
				vector3d_t N = FACE_FORWARD(sp.Ng, sp.N, wo);
				const photon_t *nearest = diffuseMap.findNearest(sp.P, N, 1.0f);
				if(nearest) col += nearest->color();
		}
		else
		{
		
		PFLOAT radius = curRadius2;    //actually the square radius... used for SPPM

		foundPhoton_t *gathered = (foundPhoton_t *)alloca(nSearch * sizeof(foundPhoton_t)); //need to be removed
		int nGathered=0;
		
		if(diffuseMap.nPhotons() > 0) nGathered = diffuseMap.gather(sp.P, gathered, nSearch, radius); //the photon gatherd function, need to be changed


		if(nGathered > 0)
		{
			if(nGathered > _nMax)
			{
				_nMax = nGathered;
				std::cout << "maximum Photons: "<<_nMax<<", radius: "<<radius<<"\n";
				if(_nMax == 10) for(int j=0; j < nGathered; ++j ) std::cout<<"col:"<<gathered[j].photon->color()<<"\n";
			}

			float scale = 1.f; //1.f / ( float(diffuseMap.nPaths()) * radius * M_PI); // scale for energy normolize
			for(int i=0; i<nGathered; ++i)
			{
				if( 1 || vector3d_t( gathered[i].photon->shadeN).normalize() * sp.Ng.normalize() > 0.5f ) // using copy constructor to removing the const effect 
				{
					curPhotons++;
					vector3d_t pdir = gathered[i].photon->direction();
					color_t surfCol = material->eval(state, sp, wo, pdir, BSDF_ALL); // need to change for BRDF
					col += surfCol * gathered[i].photon->color();// * std::fabs(sp.N*pdir); //< wrong!?
				}
			}
		}
		
		// add caustics need to be changed for sppm
		//merged into diffuseMap
		}

		state.raylevel++;
		// glossy are removed now
		if(state.raylevel <= rDepth)
		{
			//...perfect specular reflection/refraction with recursive raytracing...
			if(bsdfs & (BSDF_SPECULAR | BSDF_FILTER))
			{
				state.includeLights = true;
				bool reflect=false, refract=false;
				vector3d_t dir[2];
				color_t rcol[2], vcol;
				material->getSpecular(state, sp, wo, reflect, refract, dir, rcol);
				const volumeHandler_t *vol;
				if(reflect)
				{
					diffRay_t refRay(sp.P, dir[0], MIN_RAYDIST);
					spDiff.reflectedRay(ray, refRay); // compute the ray differentaitl
					color_t integ = color_t(integrate(state, refRay) );
					if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * refRay.dir < 0)))
					{
						if(vol->transmittance(state, refRay, vcol)) integ *= vcol;
					}
					col += integ * rcol[0];
				}
				if(refract)
				{
					diffRay_t refRay(sp.P, dir[1], MIN_RAYDIST);
					spDiff.refractedRay(ray, refRay, material->getMatIOR());
					colorA_t integ = integrate(state, refRay);
					if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * refRay.dir < 0)))
					{
						if(vol->transmittance(state, refRay, vcol)) integ *= vcol;
					}
					col += color_t(integ) * rcol[1];
					alpha = integ.A;
				}
			}			
		}
		--state.raylevel;
		
		CFLOAT m_alpha = material->getAlpha(state, sp, wo);
		alpha = m_alpha + (1.f-m_alpha)*alpha;
	}

	else //nothing hit, return background    fixed me for SPPM
	{
		if(background)
		{
			constantColor += (*background)(ray, state, false);  //this maybe wrong for sppm
		}
	}

	state.userdata = o_udat;
	state.includeLights = oldIncludeLights;
	return colorA_t(col, alpha);	
}

integrator_t* SPPM::factory(paraMap_t &params, renderEnvironment_t &render)
{
	//int shadowDepth=5; //may used 
	int raydepth=5;
	int _passNum = 1000;
	int numPhotons = 500000;
	int bounces = 5;
	float times = 1.f;

    params.getParam("raydepth", raydepth);
	params.getParam("passNums", _passNum);
	params.getParam("bounces", bounces);
	params.getParam("times", times); // initial radius times

	SPPM* ite = new SPPM(numPhotons, _passNum);
	ite->rDepth = raydepth;
	ite->nSearch = 10000;
	ite->maxBounces = bounces;
	ite->initialFactor = times;

	return ite;
}

extern "C"
{

	YAFRAYPLUGIN_EXPORT void registerPlugin(renderEnvironment_t &render)
	{
		render.registerFactory("SPPM", SPPM::factory);
	}

}
__END_YAFRAY