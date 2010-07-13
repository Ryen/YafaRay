//Quesions: 1. The Material's emit method seems called twice?
//			       2. The direct photons are not deposited now, and its strategy need double check.
//                3. The reseed strategy need more tests.

#include <integrators/sppm.h>
#include <yafraycore/scr_halton.h>
#include <sstream>
#include <cmath>
#include <algorithm>
__BEGIN_YAFRAY

const int nMaxGather = 10000; //used to gather all the photon in the radius.

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

	sDepth = 5;
	trShad = false;
	bHashgrid = false;

	//initialize the halton variable
	//hal2.setBase(2);
	//hal3.setBase(3);
	//hal5.setBase(5);
	//hal7.setBase(7);

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

	gTimer.addEvent("rendert");
	gTimer.start("rendert");
	imageFilm->init(passNum); // need to doule-check how it effect sppm.
	
	const camera_t* camera = scene->getCamera();

	maxDepth = 0.f;
	minDepth = 1e38f;

	//point3d_t a(1e38f,1e38f,1e38f);
	//point3d_t g(0.f, 0.f, 0.f); // Is it right for the initial g point?
	//bound_t bBox(a, g);

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
				//bBox.include(sp.P);
			}
		}
		if(maxDepth > 0.f) maxDepth = 1.f / (maxDepth - minDepth);
	}

	initializePPM(PM_IRE);

	renderPass(1, 0, false);

	if(PM_IRE)
	{
		totalnPhotons = 0;
		PM_IRE = false;
	}

	int size = camera->resX() * camera->resY();
	for(int i=0; i<passNum; ++i) // progress pass
	{
		if(scene->getSignals() & Y_SIG_ABORT) break;
		//imageFilm->setAAThreshold(AA_threshold);
		imageFilm->nextPass(false, integratorName);
		nRefined = 0;
		renderPass(1, 1 + (i-1)*1, false); // offset are only related to the passNum, since we alway have only one sample.
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
				HitPoint &hp = hitPoints[index];

				if(PM_IRE)
				{
					traceIRERay(rstate, c_ray, hp);
					continue;
				}

				GatherInfo gInfo = traceGatherRay(rstate, c_ray, hp);
				gInfo.photonFlux *= scene->volIntegrator->transmittance(rstate, c_ray);
				gInfo.constantRandiance *= scene->volIntegrator->transmittance(rstate, c_ray);
				gInfo.constantRandiance += scene->volIntegrator->integrate(rstate, c_ray); // Now using it to simulate for volIntegrator not using PPM

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
				if(gInfo.photonCount > 0)
				{
					const unsigned long long pcount = hp.accPhotonCount + gInfo.photonCount;
					const float g = _alpha * pcount / (hp.accPhotonCount * _alpha + gInfo.photonCount);
					hp.accPhotonCount = pcount;
					hp.accPhotonFlux = (hp.accPhotonFlux + gInfo.photonFlux) * g;
					hp.radius2 *= g;
					nRefined++;
				}

				//radiance estimate
				colorA_t color = hp.accPhotonFlux / (hp.radius2 * M_PI * totalnPhotons) + gInfo.constantRandiance;
				
				color.A = gInfo.constantRandiance.A; // This line may not needed.

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
void SPPM::prePass(int samples, int offset, bool adaptive)
{
	std::stringstream set;
	gTimer.addEvent("prePass");
	gTimer.start("prePass");

	Y_INFO << integratorName << ": Starting Photon tracing pass...\n";

	if(trShad)
	{
		set << "ShadowDepth [" << sDepth << "]";
	}
	if(!set.str().empty()) set << "+";
	set << "RayDepth [" << rDepth << "]";

	if(bHashgrid) photonGrid.clear();
	else diffuseMap.clear();


	background = scene->getBackground();
	lights = scene->lights;
	std::vector<light_t*> tmplights;

	//background do not emit photons, or it is merged into normal light?	
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
		pcol = tmplights[i]->emitPhoton(.5, .5, .5, .5, ray, lightPdf); // fix me. strange error
		lightNumPdf = lightPowerD->func[i] * lightPowerD->invIntegral;
		pcol *= fNumLights*lightPdf/lightNumPdf; //remember that lightPdf is the inverse of the pdf, hence *=...
		Y_INFO << integratorName << ": Light ["<<i+1<<"] Photon col:"<<pcol<<" | lnpdf: "<<lightNumPdf<<"\n";
	}
	
	delete[] energies;
	
	//shoot photons
	bool done=false;
	unsigned int curr=0;

	//reseed the halton sequence
	//int pass_offset = offset * 7; 
	if(offset+3 < 50)
	{    
	hal2.setBase(prims[1+offset]);
	hal3.setBase(prims[2+offset]);
	hal4.setBase(prims[3+offset]);
	}

    hal2.setStart(offset);
    hal3.setStart(offset);
    hal4.setStart(offset);

	surfacePoint_t sp;
	random_t prng(offset*(4517)+123);
	renderState_t state(&prng);
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
	//Pregather  photons

	float invDiffPhotons = 1.f / (float)nPhotons;

	unsigned int ndPhotonStored = 0;
	unsigned int ncPhotonStored = 0;
	

	while(!done)
	{
		if(scene->getSignals() & Y_SIG_ABORT) {  pb->done(); if(!intpb) delete pb; return; }
		state.chromatic = true;
		state.wavelength = scrHalton(5, curr);
		
		//s1 = RI_vdC(curr+offset);
		//s2 = scrHalton(2, curr+offset);
		//s3 = scrHalton(3, curr+offset);
		//s4 = scrHalton(4, curr+offset);
		
		//s1 = ourRandom(); //(*(state.prng))(); 
		//s2 = ourRandom(); //(*(state.prng))(); 
		//s3 = ourRandom(); //(*(state.prng))(); 
		//s4 = ourRandom(); //(*(state.prng))(); 

		//Halton
       s1 = RI_vdC(curr);
       if(offset+3 < 50)
       {
               s2 = hal2.getNext();
               s3 = hal3.getNext();
               s4 = hal4.getNext();
       }
       else
       {
               s2 = ourRandom();
               s3 = ourRandom();
               s4 = ourRandom();
       }
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
		
			if(1 || !directPhoton) // deposit all the photons 
			{
				//deposit diffuse photon on surface
				if( (!causticPhoton) && (bsdfs & (BSDF_DIFFUSE)))
				{		
					photon_t np(wi, sp.P, pcol);// pcol used here
					np.shadeN = sp.Ng;

					if(bHashgrid) photonGrid.pushPhoton(np);
					else
					{
					diffuseMap.pushPhoton(np);
					diffuseMap.setNumPaths(curr); 
					}

					ndPhotonStored++;
				}
				//deposit caustic photon on surface       This count method seems not work properly. Need to double-check
				else if(causticPhoton && (bsdfs & (BSDF_DIFFUSE | BSDF_GLOSSY)))
				{
					photon_t np(wi, sp.P, pcol);// pcol used here
					np.shadeN = sp.Ng;

					if(bHashgrid) photonGrid.pushPhoton(np);
					else
					{
						diffuseMap.pushPhoton(np);
						diffuseMap.setNumPaths(curr); 
					}

					ncPhotonStored++;
				}
			}
			//if(bsdfs & (BSDF_DIFFUSE | BSDF_GLOSSY))
			//{
			//	//deposit photon on surface
			//	if(causticPhoton)
			//	{
			//		photon_t np(wi, sp.P, pcol);
			//		causticMap.pushPhoton(np);
			//		causticMap.setNumPaths(curr);
			//		ncPhotonStored++;
			//	}
			//	else
			//	{
			//		photon_t np(wi, sp.P, pcol);
			//		diffuseMap.pushPhoton(np);
			//		diffuseMap.setNumPaths(curr);
			//		ndPhotonStored++;
			//	}
			//}
			// need to break in the middle otherwise we scatter the photon and then discard it => redundant
			if(nBounces == maxBounces)break;  

			// scatter photon
			int d5 = 3*nBounces + 5;

			s5 = scrHalton(d5, curr);
			s6 = scrHalton(d5+1, curr);
			s7 = scrHalton(d5+2, curr);

			//hal8.setBase(d5);
			//hal9.setBase(d5+1);
			//hal10.setBase(d5+2);

			//hal8.setStart(pass_offset + curr);
			//hal9.setStart(pass_offset + curr);
			//hal10.setStart(pass_offset + curr);

			//s5 =  hal8.getNext();     //ourRandom();  //(*(state.prng))(); //
			//s6 =  hal9.getNext();     //ourRandom();  //(*(state.prng))(); //
			//s7 =  hal10.getNext();   //ourRandom();  //(*(state.prng))(); //
			
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

			//// Russian Roulette
			//if(nBounces > 3)
			//{
			//	float continueProbability = 0.75f;
			//	if (ourRandom() > continueProbability) break;
			//	pcol *= 1.f / continueProbability;
			//}

		}
		++curr;
		if(curr % pbStep == 0) pb->update();
		done = (curr >= nPhotons);
	}
	pb->done();
	pb->setTag("Photon map built.");
	Y_INFO << integratorName << ":Photon map built.\n";
	Y_INFO << integratorName << ": Shot "<<curr<<" photons from " << numDLights << " light(s)\n";
	Y_INFO << integratorName << ": Stored caustic photons: " << ncPhotonStored << yendl;
	Y_INFO << integratorName << ": Stored diffuse photons: " << ndPhotonStored << yendl;
	delete lightPowerD;

	totalnPhotons +=  nPhotons;//diffuseMap.nPaths(); // used for SPPM

	Y_INFO << integratorName << ": Stored photons: "<<diffuseMap.nPhotons()<<"\n";
	
	if(bHashgrid)
	{
		Y_INFO << integratorName << ": Building photons hashgrid:\n";
		pb->setTag("Building photons hashgrid...");
		photonGrid.updateGrid();
		Y_INFO << integratorName << ": Done.\n";
	}
	else
	{
		if(diffuseMap.nPhotons() > 0) //
		{
			Y_INFO << integratorName << ": Building photons kd-tree:\n";
			pb->setTag("Building photons kd-tree...");
			diffuseMap.updateTree();
			Y_INFO << integratorName << ": Done.\n";
		}
	}

	if(diffuseMap.nPhotons() < 50) { Y_ERROR << integratorName << ": Too few photons, stopping now.\n"; return; }
	tmplights.clear();

	if(!intpb) delete pb;
	
	gTimer.stop("prePass");
	if(bHashgrid)
		Y_INFO << integratorName << ": PhotonGrid building time: " << gTimer.getTime("prePass") << "\n";
	else
		Y_INFO << integratorName << ": PhotonMap building time: " << gTimer.getTime("prePass") << "\n";

	return;
}

//collect photon information, it should return non-normalized flux for each pixel.
colorA_t SPPM::integrate(renderState_t &state, diffRay_t &ray/*, sampler_t &sam*/) const
{
	return colorA_t(0.f);	
}


GatherInfo SPPM::traceGatherRay(yafaray::renderState_t &state, yafaray::diffRay_t &ray, yafaray::HitPoint &hp)
{
	static int _nMax=0;
	static int calls=0;
	++calls;
	color_t col(0.0);
	GatherInfo gInfo;
	gInfo.constantRandiance = colorA_t(0.f);
	gInfo.photonCount = 0;
	gInfo.photonFlux = colorA_t(0.f);

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
		gInfo.constantRandiance += material->emit(state, sp, wo);
		state.includeLights = false;
		spDifferentials_t spDiff(sp, ray);

		//if(bsdfs & BSDF_DIFFUSE)
		//{
		//	gInfo.constantRandiance += estimateAllDirectLight(state, sp, wo);
		//}

		
		PFLOAT radius = hp.radius2;    //actually the square radius... used for SPPM

		foundPhoton_t *gathered = (foundPhoton_t *)alloca(nMaxGather * sizeof(foundPhoton_t)); //need to be removed
		int nGathered=0;
		
		if(bHashgrid) nGathered = photonGrid.gather(sp.P, gathered, nMaxGather, radius);
		else
		{
			if(diffuseMap.nPhotons() > 0) 
				nGathered = diffuseMap.gather(sp.P, gathered, nMaxGather, radius); //the photon gatherd function, need to be changed
		}

		if(nGathered > 0)
		{
			if(nGathered > _nMax)
			{
				_nMax = nGathered;
				std::cout << "maximum Photons: "<<_nMax<<", radius: "<<radius<<"\n";
				if(_nMax == 10) for(int j=0; j < nGathered; ++j ) std::cout<<"col:"<<gathered[j].photon->color()<<"\n";
			}
			float scale = 1.f ;// (radius * M_PI); // scale for energy normolize
			for(int i=0; i<nGathered; ++i)
			{
				if( 1 || vector3d_t( gathered[i].photon->shadeN).normalize() * sp.Ng.normalize() > 0.5f ) // using copy constructor to removing the const effect 
				{
					gInfo.photonCount++;
					vector3d_t pdir = gathered[i].photon->direction();
					color_t surfCol = material->eval(state, sp, wo, pdir, BSDF_DIFFUSE); // need to change for BRDF
					gInfo.photonFlux += surfCol * scale * gathered[i].photon->color();// * std::fabs(sp.N*pdir); //< wrong!?
				}
			}
		}
		
		// add caustics need to be changed for sppm
		//merged into diffuseMap

		state.raylevel++;
		if(state.raylevel <= rDepth)
		{
			// dispersive effects with recursive raytracing:
			if( (bsdfs & BSDF_DISPERSIVE) && state.chromatic )
			{
				state.includeLights = true; //debatable...
				int dsam = 1;  // Must trace one sample per pass or the progressive process will not work properly. Need double-check
				int oldDivision = state.rayDivision;
				int oldOffset = state.rayOffset;
				float old_dc1 = state.dc1, old_dc2 = state.dc2;
				if(state.rayDivision > 1) dsam = std::max(1, dsam/oldDivision);
				state.rayDivision *= dsam;
				int branch = state.rayDivision*oldOffset;
				float d_1 = 1.f/(float)dsam;
				float ss1 = RI_S(state.pixelSample + state.samplingOffs);
				color_t dcol(0.f), vcol(1.f);
				GatherInfo cing;
				vector3d_t wi;
				const volumeHandler_t* vol;
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
						color_t wl_col;
						wl2rgb(state.wavelength, wl_col);
						state.chromatic = false;
						refRay = diffRay_t(sp.P, wi, MIN_RAYDIST);
						cing = traceGatherRay(state, refRay, hp);
						cing.photonFlux *= mcol * wl_col;
						cing.constantRandiance *= mcol * wl_col;
						state.chromatic = true;
					}
				}
				if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * refRay.dir < 0)))
				{
					vol->transmittance(state, refRay, vcol);
					cing.photonFlux *= vcol;
					cing.constantRandiance *= vcol;
				}

				gInfo.constantRandiance += cing.constantRandiance;
				gInfo.photonFlux += cing.photonFlux;
				gInfo.photonCount += cing.photonCount;

				state.rayDivision = oldDivision;
				state.rayOffset = oldOffset;
				state.dc1 = old_dc1; state.dc2 = old_dc2;
			}
			
			// glossy reflection with recursive raytracing:  Pure GLOSSY material doesn't hold photons?
			if( bsdfs & (BSDF_GLOSSY))
			{
				state.includeLights = false;
				int gsam = 1; // Must trace one sample per pass or the progressive process will not work properly. Need double-check
				int oldDivision = state.rayDivision;
				int oldOffset = state.rayOffset;
				float old_dc1 = state.dc1, old_dc2 = state.dc2;
				if(state.rayDivision > 1) gsam = std::max(1, gsam/oldDivision);
				state.rayDivision *= gsam;
				int branch = state.rayDivision*oldOffset;
				int offs = gsam * state.pixelSample + state.samplingOffs;
				float d_1 = 1.f/(float)gsam;
				color_t gcol(0.f), vcol(1.f);
				GatherInfo ging;
				vector3d_t wi;
				const volumeHandler_t* vol;
				diffRay_t refRay;
				for(int ns=0; ns<gsam; ++ns)
				{
					state.dc1 = scrHalton(2*state.raylevel+1, branch + state.samplingOffs);
					state.dc2 = scrHalton(2*state.raylevel+2, branch + state.samplingOffs);
					state.rayOffset = branch;
					++branch;
					float s1 = RI_vdC(offs + ns);
					float s2 = scrHalton(2, offs + ns);
					if(oldDivision > 1) // create generalized halton sequence
					{
						s1 = addMod1(s1, old_dc1);
						s2 = addMod1(s2, old_dc2);
					}
					sample_t s(s1, s2, BSDF_REFLECT|BSDF_TRANSMIT|BSDF_GLOSSY);
					color_t mcol = material->sample(state, sp, wo, wi, s);
					if(s.pdf > 1.0e-5f && (s.sampledFlags & BSDF_GLOSSY))
					{
						mcol *= std::fabs(wi*sp.N)/s.pdf;
						refRay = diffRay_t(sp.P, wi, MIN_RAYDIST);
						ging += traceGatherRay(state, refRay, hp);
						ging.photonFlux *=mcol;
						ging.constantRandiance *= mcol;
					}
					
					if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * refRay.dir < 0)))
					{
						if(vol->transmittance(state, refRay, vcol)) 
						{
							ging.photonFlux *= vcol;
							ging.constantRandiance *= vcol;
						}
					}
				}
				gInfo.constantRandiance += ging.constantRandiance * d_1;
				gInfo.photonFlux += ging.photonFlux * d_1;
				gInfo.photonCount += ging.photonCount * d_1;

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
					spDiff.reflectedRay(ray, refRay); // compute the ray differentaitl
					GatherInfo refg = traceGatherRay(state, refRay, hp);
					if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * refRay.dir < 0)))
					{
						if(vol->transmittance(state, refRay, vcol)) 
						{
							refg.constantRandiance *= vcol;
							refg.photonFlux *= vcol;
						}
					}
					gInfo.constantRandiance += refg.constantRandiance * colorA_t(rcol[0]);
					gInfo.photonFlux += refg.photonFlux * colorA_t(rcol[0]);
					gInfo.photonCount += refg.photonCount;
				}
				if(refract)
				{
					diffRay_t refRay(sp.P, dir[1], MIN_RAYDIST);
					spDiff.refractedRay(ray, refRay, material->getMatIOR());
					GatherInfo refg = traceGatherRay(state, refRay, hp);
					if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * refRay.dir < 0)))
					{
						if(vol->transmittance(state, refRay, vcol)) 
						{
							refg.constantRandiance *= vcol;
							refg.photonFlux *= vcol;
						}
					}
					gInfo.constantRandiance += refg.constantRandiance * colorA_t(rcol[1]);
					gInfo.photonFlux += refg.photonFlux * colorA_t(rcol[1]);
					gInfo.photonCount += refg.photonCount;
					alpha = refg.constantRandiance.A;
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
			gInfo.constantRandiance += (*background)(ray, state, false);  //this maybe wrong for sppm
		}
	}

	state.userdata = o_udat;
	state.includeLights = oldIncludeLights;

	gInfo.constantRandiance.A = alpha;
	return gInfo;

}


void SPPM::traceIRERay(yafaray::renderState_t &state, yafaray::diffRay_t &ray, yafaray::HitPoint &hp)
{
	static int _nMax=0;
	static int calls=0;
	++calls;
	color_t col(0.0);
	GatherInfo gInfo;
	gInfo.constantRandiance = colorA_t(0.f);
	gInfo.photonCount = 0;
	gInfo.photonFlux = colorA_t(0.f);

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

		PFLOAT radius = 1.0f;    //actually the square radius... used for SPPM

		foundPhoton_t *gathered = (foundPhoton_t *)alloca(nSearch * sizeof(foundPhoton_t)); //need to be removed
		int nGathered=0;
		
		if(diffuseMap.nPhotons() > 0) nGathered = diffuseMap.gather(sp.P, gathered, nSearch, radius); //the photon gatherd function, need to be changed

		hp.radius2 = radius;
		
		// add caustics need to be changed for sppm
		//merged into diffuseMap
		}

		//state.raylevel++;
		//// glossy are removed now
		//if(state.raylevel <= rDepth)
		//{
		//	// dispersive effects with recursive raytracing:
		//	if( (bsdfs & BSDF_DISPERSIVE) && state.chromatic )
		//	{
		//		state.includeLights = true; //debatable...
		//		int dsam = 1;
		//		int oldDivision = state.rayDivision;
		//		int oldOffset = state.rayOffset;
		//		float old_dc1 = state.dc1, old_dc2 = state.dc2;
		//		if(state.rayDivision > 1) dsam = std::max(1, dsam/oldDivision);
		//		state.rayDivision *= dsam;
		//		int branch = state.rayDivision*oldOffset;
		//		float d_1 = 1.f/(float)dsam;
		//		float ss1 = RI_S(state.pixelSample + state.samplingOffs);
		//		color_t dcol(0.f), vcol(1.f);
		//		GatherInfo cing;
		//		vector3d_t wi;
		//		const volumeHandler_t* vol;
		//		diffRay_t refRay;
		//		for(int ns=0; ns<dsam; ++ns)
		//		{
		//			state.wavelength = (ns + ss1)*d_1;
		//			state.dc1 = scrHalton(2*state.raylevel+1, branch + state.samplingOffs);
		//			state.dc2 = scrHalton(2*state.raylevel+2, branch + state.samplingOffs);
		//			if(oldDivision > 1)	state.wavelength = addMod1(state.wavelength, old_dc1);
		//			state.rayOffset = branch;
		//			++branch;
		//			sample_t s(0.5f, 0.5f, BSDF_REFLECT|BSDF_TRANSMIT|BSDF_DISPERSIVE);
		//			color_t mcol = material->sample(state, sp, wo, wi, s);
		//			if(s.pdf > 1.0e-6f && (s.sampledFlags & BSDF_DISPERSIVE))
		//			{
		//				mcol *= std::fabs(wi*sp.N)/s.pdf;
		//				color_t wl_col;
		//				wl2rgb(state.wavelength, wl_col);
		//				state.chromatic = false;
		//				refRay = diffRay_t(sp.P, wi, MIN_RAYDIST);
		//				cing = traceGatherRay(state, refRay, hp);
		//				cing.photonFlux *= mcol * wl_col;
		//				cing.constantRandiance *= mcol * wl_col;
		//				state.chromatic = true;
		//			}
		//		}
		//		if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * refRay.dir < 0)))
		//		{
		//			vol->transmittance(state, refRay, vcol);
		//			cing.photonFlux *= vcol;
		//			cing.constantRandiance *= vcol;
		//		}

		//		gInfo.constantRandiance += cing.constantRandiance;
		//		gInfo.photonFlux += cing.photonFlux;
		//		gInfo.photonCount += cing.photonCount;

		//		state.rayDivision = oldDivision;
		//		state.rayOffset = oldOffset;
		//		state.dc1 = old_dc1; state.dc2 = old_dc2;
		//	}
		//	
		//	// glossy reflection with recursive raytracing:
		//	if( bsdfs & (BSDF_GLOSSY))
		//	{
		//		state.includeLights = false;
		//		int gsam = 1;
		//		int oldDivision = state.rayDivision;
		//		int oldOffset = state.rayOffset;
		//		float old_dc1 = state.dc1, old_dc2 = state.dc2;
		//		if(state.rayDivision > 1) gsam = std::max(1, gsam/oldDivision);
		//		state.rayDivision *= gsam;
		//		int branch = state.rayDivision*oldOffset;
		//		int offs = gsam * state.pixelSample + state.samplingOffs;
		//		float d_1 = 1.f/(float)gsam;
		//		color_t gcol(0.f), vcol(1.f);
		//		GatherInfo ging;
		//		vector3d_t wi;
		//		const volumeHandler_t* vol;
		//		diffRay_t refRay;
		//		for(int ns=0; ns<gsam; ++ns)
		//		{
		//			state.dc1 = scrHalton(2*state.raylevel+1, branch + state.samplingOffs);
		//			state.dc2 = scrHalton(2*state.raylevel+2, branch + state.samplingOffs);
		//			state.rayOffset = branch;
		//			++branch;
		//			float s1 = RI_vdC(offs + ns);
		//			float s2 = scrHalton(2, offs + ns);
		//			if(oldDivision > 1) // create generalized halton sequence
		//			{
		//				s1 = addMod1(s1, old_dc1);
		//				s2 = addMod1(s2, old_dc2);
		//			}
		//			sample_t s(s1, s2, BSDF_REFLECT|BSDF_TRANSMIT|BSDF_GLOSSY);
		//			color_t mcol = material->sample(state, sp, wo, wi, s);
		//			if(s.pdf > 1.0e-5f && (s.sampledFlags & BSDF_GLOSSY))
		//			{
		//				mcol *= std::fabs(wi*sp.N)/s.pdf;
		//				refRay = diffRay_t(sp.P, wi, MIN_RAYDIST);
		//				ging = traceGatherRay(state, refRay, hp);
		//				ging.photonFlux *=mcol;
		//				ging.constantRandiance *= mcol;
		//			}
		//			
		//			if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * refRay.dir < 0)))
		//			{
		//				if(vol->transmittance(state, refRay, vcol)) 
		//				{
		//					ging.photonFlux *= vcol;
		//					ging.constantRandiance *= vcol;
		//				}
		//			}
		//		}
		//		gInfo.constantRandiance += ging.constantRandiance;
		//		gInfo.photonFlux += ging.photonFlux;
		//		gInfo.photonCount += ging.photonCount;

		//		//restore renderstate
		//		state.rayDivision = oldDivision;
		//		state.rayOffset = oldOffset;
		//		state.dc1 = old_dc1; state.dc2 = old_dc2;
		//	}			
		//	//...perfect specular reflection/refraction with recursive raytracing...
		//	if(bsdfs & (BSDF_SPECULAR | BSDF_FILTER))
		//	{
		//		state.includeLights = true;
		//		bool reflect=false, refract=false;
		//		vector3d_t dir[2];
		//		color_t rcol[2], vcol;
		//		material->getSpecular(state, sp, wo, reflect, refract, dir, rcol);
		//		const volumeHandler_t *vol;
		//		if(reflect)
		//		{
		//			diffRay_t refRay(sp.P, dir[0], MIN_RAYDIST);
		//			spDiff.reflectedRay(ray, refRay); // compute the ray differentaitl
		//			GatherInfo refg = traceGatherRay(state, refRay, hp);
		//			if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * refRay.dir < 0)))
		//			{
		//				if(vol->transmittance(state, refRay, vcol)) 
		//				{
		//					refg.constantRandiance *= vcol;
		//					refg.photonFlux *= vcol;
		//				}
		//			}
		//			gInfo.constantRandiance += refg.constantRandiance * colorA_t(rcol[0]);
		//			gInfo.photonFlux += refg.photonFlux * colorA_t(rcol[0]);
		//			gInfo.photonCount += refg.photonCount;
		//		}
		//		if(refract)
		//		{
		//			diffRay_t refRay(sp.P, dir[1], MIN_RAYDIST);
		//			spDiff.refractedRay(ray, refRay, material->getMatIOR());
		//			GatherInfo refg = traceGatherRay(state, refRay, hp);
		//			if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * refRay.dir < 0)))
		//			{
		//				if(vol->transmittance(state, refRay, vcol)) 
		//				{
		//					refg.constantRandiance *= vcol;
		//					refg.photonFlux *= vcol;
		//				}
		//			}
		//			gInfo.constantRandiance += refg.constantRandiance * colorA_t(rcol[1]);
		//			gInfo.photonFlux += refg.photonFlux * colorA_t(rcol[1]);
		//			gInfo.photonCount += refg.photonCount;
		//			alpha = refg.constantRandiance.A;
		//		}
		//	}			
		//}
		//--state.raylevel;
		
		//CFLOAT m_alpha = material->getAlpha(state, sp, wo);
		//alpha = m_alpha + (1.f-m_alpha)*alpha;
	}

	else //nothing hit, return background    fixed me for SPPM
	{
		if(background)
		{
			//gInfo.constantRandiance += (*background)(ray, state, false);  //this maybe wrong for sppm
		}
	}

	state.userdata = o_udat;
	state.includeLights = oldIncludeLights;

	return;	
}


void SPPM::initializePPM(bool us_PM)
{
	const camera_t* camera = scene->getCamera();
	hitPoints.reserve(camera->resX() * camera->resY()); // used for SPPM
	bound_t bBox = scene->getSceneBound(); // Now using Scene Bound 

	// initialize SPPM statistics
	float initialRadius = ((bBox.longX() + bBox.longY() + bBox.longZ()) / 3.f) / ((camera->resX() + camera->resY()) / 2.0f) * 2.f ;

	initialRadius = std::min(initialRadius, 1.f); //Fix the overflow bug

	int size = camera->resX() * camera->resY();
	for(int i = 0; i < size; i++)
	{
		HitPoint ts;

		ts.accPhotonFlux  = colorA_t(0.f);
		ts.accPhotonCount = 0;
		if(!PM_IRE) ts.radius2 = (initialRadius * initialFactor) * (initialRadius * initialFactor);

		ts.acclightsource = colorA_t(0.f);
		ts.constanthits = 0;
		ts.surfacehits = 0;

		hitPoints.push_back(ts);
	}
	
	if(bHashgrid) photonGrid.setParm(initialRadius*2.f, nPhotons, bBox);

}

integrator_t* SPPM::factory(paraMap_t &params, renderEnvironment_t &render)
{
	bool transpShad=false;
	int shadowDepth=5; //may used when integrate Direct Light
	int raydepth=5;
	int _passNum = 1000;
	int numPhotons = 500000;
	int bounces = 5;
	float times = 1.f;

	params.getParam("transpShad", transpShad);
	params.getParam("shadowDepth", shadowDepth);
	params.getParam("raydepth", raydepth);
	params.getParam("photons", numPhotons);
	params.getParam("passNums", _passNum);
	params.getParam("bounces", bounces);
	params.getParam("times", times); // initial radius times

	SPPM* ite = new SPPM(numPhotons, _passNum);
	ite->sDepth = shadowDepth;
	ite->trShad = transpShad;
	ite->rDepth = raydepth;
	ite->nSearch = 100;
	ite->maxBounces = bounces;
	ite->initialFactor = times;
	ite->PM_IRE = false;

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