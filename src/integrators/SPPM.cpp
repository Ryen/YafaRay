#include "SPPM.h"
#include <sstream>
#include <cmath>
#include <algorithm>
__BEGIN_YAFRAY


SPPM::SPPM(unsigned int dPhotons, unsigned int cPhotons, int _passnum)
{
	nPhotons = dPhotons;
	nSearch = dPhotons; // need to find all photons in the region.
	nCausPhotons = cPhotons;
	nCausSearch = cPhotons; // same as above
	passNum = _passnum;
	totalnPhotons = 0;
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


	if(intpb) intpb->setTag(passString.str().c_str());

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

	// initialize SPPM statistics
	float initialRadius = ((bBox.longX() + bBox.longY() + bBox.longZ()) / 3.f) / ((camera->resX() + camera->resY()) / 2.0f) * 2.f ;
	
	std::vector<shared_statistics>::iterator itr;
	for(itr = progressiveData.begin(); itr != progressiveData.end(); itr++)
	{
		itr->accflux = colorA_t(0.f);
		itr->photoncount = 0;
		itr->radius = initialRadius;
	}
	
  //  //
	//renderPass(AA_samples, 0, false);

	for(int i=0; i<passNum; ++i) // progress pass
	{
		if(scene->getSignals() & Y_SIG_ABORT) break;
		//imageFilm->setAAThreshold(AA_threshold);
		imageFilm->nextPass(false, integratorName);
		renderPass(1, 1 + (i-1)*1, false); // offset seems only used for sampling?
	}
	maxDepth = 0.f;
	gTimer.stop("rendert");
	Y_INFO << integratorName << ": Overall rendertime: "<< gTimer.getTime("rendert")<<"s\n";

	return true;
}

bool SPPM::renderPass(int samples,int offset, bool adaptive)
{
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

			if(adaptive) // now not need for sppm
			{
				if(!imageFilm->doMoreSamples(j, i)) continue;
			}

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
				if(passNum>1) // ALWAYS True for SPPM
				{
					dx = RI_S(rstate.pixelSample, rstate.samplingOffs);
					dy = RI_vdC(rstate.pixelSample, rstate.samplingOffs);
				}
				else if(n_samples > 1) //not need now
				{
					dx = (0.5+(PFLOAT)sample)*d1;
					dy = RI_LP(sample+rstate.samplingOffs);
				}
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

				//for sppm progressive
				int index = i*camera->resX() + j; //not possible, w is a protected member.
				shared_statistics &temp = progressiveData[index];

				curPhotons = 0;
				curRadius = temp.radius;
				colorA_t flux = integrate(rstate, c_ray); // L_o
				flux *= scene->volIntegrator->transmittance(rstate, c_ray); // T
				flux += scene->volIntegrator->integrate(rstate, c_ray); // L_v

				// progressive refinement
				const float _alpha = 0.8;
				float g = std::min((temp.photoncount + _alpha * curPhotons) / (temp.photoncount + curPhotons), 1.0f);
				temp.radius = temp.radius * sqrt(g);
				temp.photoncount += curPhotons * _alpha;
				temp.accflux = (temp.accflux + flux) * g;				
				

				//radiance estimate
				colorA_t color = temp.accflux / (temp.radius * temp.radius * 3.141592 * totalnPhotons);				
				imageFilm->addSample(wt * color, j, i, dx, dy, &a);
				
				
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
bool SPPM::prepass(int samples, int offset, bool adaptive)
{
	std::stringstream set;
	gTimer.addEvent("prepass");
	gTimer.start("prepass");

	Y_INFO << integratorName << ": Starting preprocess...\n";

	//if(trShad)
	//{
	//	set << "ShadowDepth [" << sDepth << "]";
	//}
	//if(!set.str().empty()) set << "+";
	//set << "RayDepth [" << rDepth << "]";

	diffuseMap.clear();
	causticMap.clear();
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
		if(lights[i]->shootsDiffuseP())
		{
			numDLights++;
			tmplights.push_back(lights[i]);
		}
	}
	
	fNumLights = (float)numDLights;
	energies = new float[numDLights];

	for(int i=0;i<numDLights;++i) energies[i] = tmplights[i]->totalEnergy().energy();

	lightPowerD = new pdf1D_t(energies, numDLights);
	
	Y_INFO << integratorName << ": Light(s) photon color testing for diffuse map:\n";
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
	if(intpb) pb = intpb;
	else pb = new ConsoleProgressBar_t(80);
	
	Y_INFO << integratorName << ": Building diffuse photon map...\n";
	
	pb->init(128);
	pbStep = std::max(1U, nPhotons/128);
	pb->setTag("Building diffuse photon map...");
	//Pregather diffuse photons

	float invDiffPhotons = 1.f / (float)nPhotons;
	
	while(!done)
	{
		if(scene->getSignals() & Y_SIG_ABORT) {  pb->done(); if(!intpb) delete pb; return false; }
		//state.chromatic = true;
		//state.wavelength = scrHalton(5, curr);

		s1 = RI_vdC(curr);
		s2 = scrHalton(2, curr);
		s3 = scrHalton(3, curr);
		s4 = scrHalton(4, curr);

		sL = float(curr) * invDiffPhotons;
		int lightNum = lightPowerD->DSample(sL, &lightNumPdf);
		if(lightNum >= numDLights){ Y_ERROR << integratorName << ": lightPDF sample error! "<<sL<<"/"<<lightNum<<"... stopping now.\n"; delete lightPowerD; return false; }

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
				if(!causticPhoton)
				{
					photon_t np(wi, sp.P, pcol);// pcol used here
					diffuseMap.pushPhoton(np);
					diffuseMap.setNumPaths(curr);
				}

			}
			// need to break in the middle otherwise we scatter the photon and then discard it => redundant
			if(nBounces == maxBounces) break;
			// scatter photon
			int d5 = 3*nBounces + 5;

			s5 = scrHalton(d5, curr);
			s6 = scrHalton(d5+1, curr);
			s7 = scrHalton(d5+2, curr);
			
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
		}
		++curr;
		if(curr % pbStep == 0) pb->update();
		done = (curr >= nPhotons);
	}
	pb->done();
	pb->setTag("Diffuse photon map built.");
	Y_INFO << integratorName << ": Diffuse photon map built.\n";
	Y_INFO << integratorName << ": Shot "<<curr<<" photons from " << numDLights << " light(s)\n";

	delete lightPowerD;

	tmplights.clear();

	for(int i=0;i<(int)lights.size();++i)
	{
		if(lights[i]->shootsCausticP())
		{
			numCLights++;
			tmplights.push_back(lights[i]);
		}
	}

	if(numCLights > 0)
	{
		
		done = false;
		curr=0;

		fNumLights = (float)numCLights;
		energies = new float[numCLights];

		for(int i=0;i<numCLights;++i) energies[i] = tmplights[i]->totalEnergy().energy();

		lightPowerD = new pdf1D_t(energies, numCLights);
		
		Y_INFO << integratorName << ": Light(s) photon color testing for caustics map:\n";
		for(int i=0;i<numCLights;++i)
		{
			pcol = tmplights[i]->emitPhoton(.5, .5, .5, .5, ray, lightPdf);
			lightNumPdf = lightPowerD->func[i] * lightPowerD->invIntegral;
			pcol *= fNumLights*lightPdf/lightNumPdf; //remember that lightPdf is the inverse of the pdf, hence *=...
			Y_INFO << integratorName << ": Light ["<<i+1<<"] Photon col:"<<pcol<<" | lnpdf: "<<lightNumPdf<<"\n";
		}
		totalnPhotons +=  nPhotons; // used for SPPM

		delete[] energies;

		Y_INFO << integratorName << ": Building caustics photon map...\n";
		pb->init(128);
		pbStep = std::max(1U, nCausPhotons / 128);
		pb->setTag("Building caustics photon map...");
		//Pregather caustic photons
		
		float invCaustPhotons = 1.f / (float)nCausPhotons;
		
		while(!done)
		{
			if(scene->getSignals() & Y_SIG_ABORT) { pb->done(); if(!intpb) delete pb; return false; }
			state.chromatic = true;
			state.wavelength = scrHalton(5,curr);

			s1 = RI_vdC(curr);
			s2 = scrHalton(2, curr);
			s3 = scrHalton(3, curr);
			s4 = scrHalton(4, curr);

			sL = float(curr) * invCaustPhotons;
			int lightNum = lightPowerD->DSample(sL, &lightNumPdf);
			if(lightNum >= numCLights){ Y_ERROR << integratorName << ": lightPDF sample error! "<<sL<<"/"<<lightNum<<"... stopping now.\n"; delete lightPowerD; return false; }

			pcol = tmplights[lightNum]->emitPhoton(s1, s2, s3, s4, ray, lightPdf);
			ray.tmin = MIN_RAYDIST;
			ray.tmax = -1.0;
			pcol *= fNumLights*lightPdf/lightNumPdf; //remember that lightPdf is the inverse of th pdf, hence *=...
			if(pcol.isBlack())
			{
				++curr;
				done = (curr >= nCausPhotons);
				continue;
			}
			int nBounces=0;
			bool causticPhoton = false;
			bool directPhoton = true;
			const material_t *material = NULL;
			BSDF_t bsdfs;

			while( scene->intersect(ray, sp) )
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

				if(bsdfs & (BSDF_DIFFUSE | BSDF_GLOSSY))
				{
					if(causticPhoton)
					{
						photon_t np(wi, sp.P, pcol);
						causticMap.pushPhoton(np);
						causticMap.setNumPaths(curr);
					}
				}
				// need to break in the middle otherwise we scatter the photon and then discard it => redundant
				if(nBounces == maxBounces) break;
				// scatter photon
				int d5 = 3*nBounces + 5;

				s5 = scrHalton(d5, curr);
				s6 = scrHalton(d5+1, curr);
				s7 = scrHalton(d5+2, curr);

				pSample_t sample(s5, s6, s7, BSDF_ALL, pcol, transm);

				bool scattered = material->scatterPhoton(state, sp, wi, wo, sample);
				if(!scattered) break; //photon was absorped.

				pcol = sample.color;

				causticPhoton = ((sample.sampledFlags & (BSDF_GLOSSY | BSDF_SPECULAR | BSDF_DISPERSIVE)) && directPhoton) ||
								((sample.sampledFlags & (BSDF_GLOSSY | BSDF_SPECULAR | BSDF_FILTER | BSDF_DISPERSIVE)) && causticPhoton);
				directPhoton = (sample.sampledFlags & BSDF_FILTER) && directPhoton;
				if(state.chromatic && (sample.sampledFlags & BSDF_DISPERSIVE))
				{
					state.chromatic=false;
					color_t wl_col;
					wl2rgb(state.wavelength, wl_col);
					pcol *= wl_col;
				}
				ray.from = sp.P;
				ray.dir = wo;
				ray.tmin = MIN_RAYDIST;
				ray.tmax = -1.0;
				++nBounces;
			}
			++curr;
			if(curr % pbStep == 0) pb->update();
			done = (curr >= nCausPhotons);
		}
		
		pb->done();
		pb->setTag("Caustics photon map built.");
		delete lightPowerD;
	}
	else
	{
		Y_INFO << integratorName << ": No caustic source lights found, skiping caustic gathering...\n";		
	}
	
	Y_INFO << integratorName << ": Shot "<<curr<<" caustic photons from " << numCLights <<" light(s).\n";
	Y_INFO << integratorName << ": Stored caustic photons: "<<causticMap.nPhotons()<<"\n";
	Y_INFO << integratorName << ": Stored diffuse photons: "<<diffuseMap.nPhotons()<<"\n";
	
	if(diffuseMap.nPhotons() > 0) //
	{
		Y_INFO << integratorName << ": Building diffuse photons kd-tree:\n";
		pb->setTag("Building diffuse photons kd-tree...");
		diffuseMap.updateTree();
		Y_INFO << integratorName << ": Done.\n";
	}

	if(causticMap.nPhotons() > 0)
	{
		Y_INFO << integratorName << ": Building caustic photons kd-tree:\n";
		pb->setTag("Building caustic photons kd-tree...");
		causticMap.updateTree();
		Y_INFO << integratorName << ": Done.\n";
	}

	if(diffuseMap.nPhotons() < 50) { Y_ERROR << integratorName << ": Too few diffuse photons, stopping now.\n"; return false; }
	
	
	tmplights.clear();

	if(!intpb) delete pb;
	

	gTimer.stop("prepass");
	Y_INFO << integratorName << ": Photonmap building time: " << gTimer.getTime("prepass") << "\n";

	return true;
}



//eye pass, it should return non-normalized flux for each pixel.
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
		col += material->emit(state, sp, wo);
		state.includeLights = false;
		spDifferentials_t spDiff(sp, ray);

        // remove FG here
		foundPhoton_t *gathered = (foundPhoton_t *)alloca(nSearch * sizeof(foundPhoton_t)); //need to be removed
		PFLOAT radius = curRadius;    //actually the square radius... used for SPPM

		int nGathered=0;
		
		if(diffuseMap.nPhotons() > 0) nGathered = diffuseMap.gather(sp.P, gathered, nSearch, radius); //the photon gatherd function, need to be changed
		color_t sum(0.0);
		if(nGathered > 0)
		{
			if(nGathered > _nMax)
			{
				_nMax = nGathered;
				std::cout << "maximum Photons: "<<_nMax<<", radius: "<<radius<<"\n";
				if(_nMax == 10) for(int j=0; j < nGathered; ++j ) std::cout<<"col:"<<gathered[j].photon->color()<<"\n";
			}

			float scale = 1.f; //1.f / ( float(diffuseMap.nPaths()) * radius * M_PI); // scale for energy normolize, need to be changed!     see pbrt
			for(int i=0; i<nGathered; ++i)
			{
				vector3d_t pdir = gathered[i].photon->direction();
				color_t surfCol = material->eval(state, sp, wo, pdir, BSDF_DIFFUSE);
				col += surfCol * scale * gathered[i].photon->color();// * std::fabs(sp.N*pdir); //< wrong!?
			}
		}

		// add caustics  need to be changed for sppm
		if(bsdfs & (BSDF_DIFFUSE))
		{
			//The linker could not find this method now, why?
			//col += estimatePhotons(state, sp, causticMap, wo, nCausSearch, cRadius); // not need to scale energy, need to be changed here? 
		}
		
		state.raylevel++;
		
		// need to be checked and changed for sppm
		if(state.raylevel <= rDepth)
		{
			// dispersive effects with recursive raytracing:
			if( (bsdfs & BSDF_DISPERSIVE) && state.chromatic)
			{
				state.includeLights = false; //debatable...
				int dsam = 8;
				int oldDivision = state.rayDivision;
				int oldOffset = state.rayOffset;
				float old_dc1 = state.dc1, old_dc2 = state.dc2;
				if(state.rayDivision > 1) dsam = std::max(1, dsam/oldDivision);
				state.rayDivision *= dsam;
				int branch = state.rayDivision*oldOffset;
				float d_1 = 1.f/(float)dsam;
				float ss1 = RI_S(state.pixelSample + state.samplingOffs);
				color_t dcol(0.f), vcol(1.f);
				vector3d_t wi;
				const volumeHandler_t *vol;
				diffRay_t refRay;
				for(int ns=0; ns<dsam; ++ns)
				{
					state.wavelength = (ns + ss1)*d_1;
					state.dc1 = scrHalton(2*state.raylevel+1, branch + state.samplingOffs);
					state.dc2 = scrHalton(2*state.raylevel+2, branch + state.samplingOffs);
					if(oldDivision > 1)	state.wavelength = addMod1(state.wavelength, old_dc1);
					state.rayOffset = branch;
					++branch;
					sample_t s(0.5f, 0.5f, BSDF_REFLECT|BSDF_TRANSMIT|BSDF_DISPERSIVE);
					color_t mcol = material->sample(state, sp, wo, wi, s);
					if(s.pdf > 1.0e-6f && (s.sampledFlags & BSDF_DISPERSIVE))
					{
						mcol *= std::fabs(wi*sp.N)/s.pdf;
						state.chromatic = false;
						color_t wl_col;
						wl2rgb(state.wavelength, wl_col);
						refRay = diffRay_t(sp.P, wi, MIN_RAYDIST);
						dcol += (color_t)integrate(state, refRay) * mcol * wl_col;
						state.chromatic = true;
					}
				}
				if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.N * refRay.dir < 0)))
				{
					vol->transmittance(state, refRay, vcol);
					dcol *= vcol;
				}
				col += dcol * d_1;
				state.rayDivision = oldDivision;
				state.rayOffset = oldOffset;
				state.dc1 = old_dc1; state.dc2 = old_dc2;
			}
			// glossy reflection with recursive raytracing:
			if( bsdfs & BSDF_GLOSSY )
			{
				state.includeLights = false;
				int gsam = 8;
				int oldDivision = state.rayDivision;
				int oldOffset = state.rayOffset;
				float old_dc1 = state.dc1, old_dc2 = state.dc2;
				if(state.rayDivision > 1) gsam = std::max(1, gsam/oldDivision);
				state.rayDivision *= gsam;
				int branch = state.rayDivision*oldOffset;
				unsigned int offs = gsam * state.pixelSample + state.samplingOffs;
				float d_1 = 1.f/(float)gsam;
				color_t gcol(0.f), vcol(1.f);
				vector3d_t wi;
				const volumeHandler_t *vol;
				diffRay_t refRay;
				for(int ns=0; ns<gsam; ++ns) // This need to double check for its rightness for sppm
				{
					state.dc1 = scrHalton(2*state.raylevel+1, branch + state.samplingOffs);
					state.dc2 = scrHalton(2*state.raylevel+2, branch + state.samplingOffs);
					state.rayOffset = branch;
					++offs;
					++branch;
					
					float s1 = RI_vdC(offs);
					float s2 = scrHalton(2, offs);
					
					if(oldDivision > 1) // create generalized halton sequence
					{
						s1 = addMod1(s1, old_dc1);
						s2 = addMod1(s2, old_dc2);
					}
					
					sample_t s(s1, s2, BSDF_REFLECT|BSDF_TRANSMIT|BSDF_GLOSSY);
					color_t mcol = material->sample(state, sp, wo, wi, s);
					if(s.pdf > 1.0e-6f && (s.sampledFlags & BSDF_GLOSSY))
					{
						mcol *= std::fabs(wi*sp.N)/s.pdf;
						refRay = diffRay_t(sp.P, wi, MIN_RAYDIST);
						gcol += (color_t)integrate(state, refRay) * mcol; // need to double check it suits for sppm
					}
					if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.N * refRay.dir < 0)))
					{
						if(vol->transmittance(state, refRay, vcol)) gcol *= vcol;
					}
				}

				col += gcol * d_1; // d_1 is a normalize term. 
				//restore renderstate
				state.rayDivision = oldDivision;
				state.rayOffset = oldOffset;
				state.dc1 = old_dc1; state.dc2 = old_dc2;
			}
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
					spDiff.reflectedRay(ray, refRay);
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
	else //nothing hit, return background   
	{
		if(background)
		{
			col += (*background)(ray, state, false);  //this maybe wrong for sppm
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
	int _passNum = 10000;
	int numPhotons = 100000;
	int numCPhotons = 500000;
	int search = 100000; //need consider carefully
	int caustic_mix = 500000; // same as above
	int bounces = 5;

	//float cuRad = 23


	params.getParam("raydepth", raydepth);
	params.getParam("photons", numPhotons);
	params.getParam("cPhotons", numCPhotons);
	params.getParam("nPassNum", _passNum);
	//params.getParam("diffuseRadius", dsRad);
	//params.getParam("causticRadius", cRad);
	//params.getParam("radius", curRad);
	//params.getParam("search", search); // FIXED
	//caustic_mix = search;
	params.getParam("caustic_mix", caustic_mix);
	params.getParam("bounces", bounces);


	
	SPPM* ite = new SPPM(numPhotons, numCPhotons, _passNum);
	ite->rDepth = raydepth;
	ite->nSearch = search;
	//ite->curRadius = curRad;
	ite->nCausSearch = caustic_mix;
	ite->maxBounces = bounces;

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