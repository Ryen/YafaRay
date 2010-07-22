 //Quesions:  1.  The tracegatherRay method's russian roulette need double-check. (It seems not be theoretical right)
//                  2. Need  a correct reseed strategy to make sppm more efficient.

#include <integrators/sppm.h>
#include <yafraycore/scr_halton.h>
#include <sstream>
#include <cmath>
#include <algorithm>
__BEGIN_YAFRAY

const int nMaxGather = 10000; //used to gather all the photon in the radius. seems could get a better way to do that
SPPM::SPPM(unsigned int dPhotons, int _passnum, bool transpShad, int shadowDepth)
{
	type = SURFACE;
	intpb = 0;
	integratorName = "SPPM";
	integratorShortName = "SPPM";
	nPhotons = dPhotons;
	passNum = _passnum;
	totalnPhotons = 0;
	initialFactor = 1.f;

	sDepth = shadowDepth;
	trShad = transpShad;
	bHashgrid = false;
	
	hal1.setBase(2);
	hal2.setBase(3);
	hal3.setBase(5);
	hal4.setBase(7);

	hal1.setStart(0);
	hal2.setStart(0);
	hal3.setStart(0);
	hal4.setStart(0);
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

	if(scene->doDepth()) precalcDepths();

	initializePPM(PM_IRE); // seems could integrate into the preRender

	passString << "Rendering pass "<< 1<<" of " << std::max(1, passNum) << "...";
	if(intpb) intpb->setTag(passString.str().c_str());
	renderPass(1, 0, false);
	PM_IRE = false;

	int hpNum = camera->resX() * camera->resY();
	for(int i=0; i<passNum; ++i) //progress pass
	{
		if(scene->getSignals() & Y_SIG_ABORT) break;
		imageFilm->nextPass(false, integratorName);
		nRefined = 0;

		passString << "Rendering pass "<< i+1<<" of " << std::max(1, passNum) << "...";
		if(intpb) intpb->setTag(passString.str().c_str());
		renderPass(1, 1 + (i-1)*1, false); // offset are only related to the passNum, since we alway have only one sample.
		Y_INFO<<  integratorName <<": This pass refined "<<nRefined<<" of "<<hpNum<<" pixels."<<"\n";
	}
	maxDepth = 0.f;
	gTimer.stop("rendert");
	Y_INFO << integratorName << ": Overall rendertime: "<< gTimer.getTime("rendert")<<"s\n";
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

				GatherInfo gInfo = traceGatherRay(rstate, c_ray, hp);
				gInfo.photonFlux *= scene->volIntegrator->transmittance(rstate, c_ray);
				gInfo.constantRandiance *= scene->volIntegrator->transmittance(rstate, c_ray);
				gInfo.constantRandiance += scene->volIntegrator->integrate(rstate, c_ray); // Now using it to simulate for volIntegrator not using PPM, need more tests

				// progressive refinement
				const float _alpha = 0.7; // change to 0.8 to do a tests
				float g = 1.0f;

				// The author's refine formular
				if(gInfo.photonCount > 0)
				{
					float g = std::min((hp.accPhotonCount + _alpha * gInfo.photonCount) / (hp.accPhotonCount+ gInfo.photonCount), 1.0f);
					hp.radius2 *= g;
					hp.accPhotonCount += gInfo.photonCount * _alpha;
					hp.accPhotonFlux = (hp.accPhotonFlux + gInfo.photonFlux) * g;
					nRefined++;
				}

				//radiance estimate
				colorA_t color = hp.accPhotonFlux / (hp.radius2 * M_PI * totalnPhotons) + gInfo.constantRandiance;
				color.A = gInfo.constantRandiance.A; // maintain the alpha value, need more tests.

				imageFilm->addSample(wt * color, j, i, dx, dy, &a);

				if(do_depth)
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
		pcol = tmplights[i]->emitPhoton(.5, .5, .5, .5, ray, lightPdf); 
		lightNumPdf = lightPowerD->func[i] * lightPowerD->invIntegral;
		pcol *= fNumLights*lightPdf/lightNumPdf; //remember that lightPdf is the inverse of the pdf, hence *=...
		Y_INFO << integratorName << ": Light ["<<i+1<<"] Photon col:"<<pcol<<" | lnpdf: "<<lightNumPdf<<"\n";
	}
	
	delete[] energies;
	
	//shoot photons
	bool done=false;
	unsigned int curr=0;

	////reseed
	//if(offset+3 < 50)
	//{    
	//hal2.setBase(prims[1+offset]);
	//hal3.setBase(prims[2+offset]);
	//hal4.setBase(prims[3+offset]);
	//}

 //   hal2.setStart(offset);
 //   hal3.setStart(offset);
 //   hal4.setStart(offset);


	surfacePoint_t sp;
	random_t prng(offset*(4517)+123);
	renderState_t state(&prng);
	unsigned char userdata[USER_DATA_SIZE+7];
	state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
	state.cam = scene->getCamera();
	progressBar_t *pb;
	int pbStep;
	if(intpb) pb = intpb;
	else pb = new ConsoleProgressBar_t(80);
	
	if(bHashgrid) Y_INFO << integratorName << ": Building photon hashgrid...\n";
	else Y_INFO << integratorName << ": Building photon map...\n";
	
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

       s1 = ourRandom();//hal1.getNext();
	   s2 = ourRandom();//hal2.getNext();
       s3 = ourRandom();//hal3.getNext();
       s4 = ourRandom();//hal4.getNext();

	   //using Halton for first 50 pass, then using random(), I think a here need a better way to do it
       //s1 = RI_vdC(curr);
       //if(offset+3 < 50)
       //{
       //        s2 = hal2.getNext();
       //        s3 = hal3.getNext();
       //        s4 = hal4.getNext();
       //}
       //else
       //{
       //        s2 = ourRandom();
       //        s3 = ourRandom();
       //        s4 = ourRandom();
       //}
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
		

			//deposit photon on diffuse surface, now we only have one map for all
			if(bsdfs & (BSDF_DIFFUSE))
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
			
			// need to break in the middle otherwise we scatter the photon and then discard it => redundant
			if(nBounces == maxBounces)break;  

			// scatter photon
			int d5 = 3*nBounces + 5;

			s5 = scrHalton(d5, curr);
			s6 = scrHalton(d5+1, curr);
			s7 = scrHalton(d5+2, curr);

			pSample_t sample(s5, s6, s7, BSDF_ALL, pcol, transm);

			bool scattered = material->scatterPhoton(state, sp, wi, wo, sample);
			if(!scattered) break; //photon was absorped.  actually based on russian roulette

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
	pb->setTag("Photon map built.");
	Y_INFO << integratorName << ":Photon map built.\n";
	Y_INFO << integratorName << ": Shot "<<curr<<" photons from " << numDLights << " light(s)\n";
	Y_INFO << integratorName << ": Stored photon vertexs: " << ndPhotonStored << yendl;
	delete lightPowerD;

	totalnPhotons +=  nPhotons;	// accumulate the total photon number, not using nPath for the case of hashgrid.

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
		if(diffuseMap.nPhotons() < 50) { Y_ERROR << integratorName << ": Too few photons, stopping now.\n"; return; }
	}

	tmplights.clear();

	if(!intpb) delete pb;
	
	gTimer.stop("prePass");

	if(bHashgrid)
		Y_INFO << integratorName << ": PhotonGrid building time: " << gTimer.getTime("prePass") << "\n";
	else
		Y_INFO << integratorName << ": PhotonMap building time: " << gTimer.getTime("prePass") << "\n";

	return;
}

//now its a dummy function
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
		
		// use russian roulette to only choose one bsdf to sample
		BSDF_t  choosen_bsdf;

		//check how many
		BSDF_t collections[4];
		unsigned bcount = 0;

		if( bsdfs & BSDF_SPECULAR | BSDF_FILTER) collections[bcount++] = BSDF_SPECULAR | BSDF_FILTER;
		if( bsdfs & BSDF_GLOSSY)	collections[bcount++] = BSDF_GLOSSY;
		if( bsdfs & BSDF_DIFFUSE)	collections[bcount++] = BSDF_DIFFUSE;
		if( bsdfs & BSDF_DISPERSIVE)	collections[bcount++] = BSDF_DISPERSIVE;

		
		float p = ourRandom();

		for(int i = 0; i < bcount; i++)
		{
			if(p < (i+1)/ float(bcount))
			{
				choosen_bsdf = collections[i];
				break;
			}
		}
		bsdfs = (bsdfs & BSDF_VOLUMETRIC) ? BSDF_VOLUMETRIC : 0;
		bsdfs |= choosen_bsdf;
		
		
		PFLOAT radius = hp.radius2;    //actually the square radius... used for SPPM

		foundPhoton_t *gathered = (foundPhoton_t *)alloca(nMaxGather * sizeof(foundPhoton_t)); 

		int nGathered=0;
		
		if(bHashgrid) nGathered = photonGrid.gather(sp.P, gathered, nMaxGather, radius);
		else
		{
			if(diffuseMap.nPhotons() > 0)
			{
				if(PM_IRE)
				{
					radius = 1.0f;
					nGathered = diffuseMap.gather(sp.P, gathered, nSearch, radius);
					hp.radius2 = radius;
				}
				else
					nGathered = diffuseMap.gather(sp.P, gathered, nMaxGather, radius); //we always collected all the photon inside the radius
			}
		}

		if(nGathered > 0)
		{
			if(nGathered > _nMax)
			{
				_nMax = nGathered;
				std::cout << "maximum Photons: "<<_nMax<<", radius: "<<radius<<"\n";
				if(_nMax == 10) for(int j=0; j < nGathered; ++j ) std::cout<<"col:"<<gathered[j].photon->color()<<"\n";
			}
			//float scale = 1.f; // scale is useless now
			for(int i=0; i<nGathered; ++i)
			{
				gInfo.photonCount++;
				vector3d_t pdir = gathered[i].photon->direction();
				color_t surfCol = material->eval(state, sp, wo, pdir, BSDF_DIFFUSE); // seems could speed up using rho, (something pbrt made)
				gInfo.photonFlux += surfCol * gathered[i].photon->color();// * std::fabs(sp.N*pdir); //< wrong!?
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
				int dsam = 1;  // Must trace only one sample per pass or the progressive process will not work properly. Need double-check
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
				int gsam = 1; // Must trace only one sample per pass or the progressive process will not work properly. Need double-check
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

				//if(reflect && refract)
				//{
				//	float p = ourRandom(); // also russian roulette
				//	if(p <= 0.5) refract = false;
				//	else reflect = false;
				//}
			
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

	else //nothing hit, return background
	{
		if(background)
		{
			gInfo.constantRandiance += (*background)(ray, state, false); 
		}
	}

	state.userdata = o_udat;
	state.includeLights = oldIncludeLights;

	gInfo.constantRandiance.A = alpha; // a small trick for just hold the alpha value.

	//gInfo.photonFlux *= 4;
	//gInfo.constantRandiance *= 4;
	//gInfo.photonCount *= 4;
	return gInfo;
}

// Now not integrated
void SPPM::traceIRERay(yafaray::renderState_t &state, yafaray::diffRay_t &ray, yafaray::HitPoint &hp)
{
	static int _nMax=0;
	static int calls=0;
	++calls;

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
		vector3d_t wo = -ray.dir;
		const material_t *material = sp.material;
		material->initBSDF(state, sp, bsdfs);
		state.includeLights = false;
		spDifferentials_t spDiff(sp, ray);

		if(bsdfs & BSDF_DIFFUSE)
		{
			PFLOAT radius = 1.0f;   // need to be nsearch,
			foundPhoton_t *gathered = (foundPhoton_t *)alloca(nSearch * sizeof(foundPhoton_t)); 
			int nGathered=0;
			
			if(diffuseMap.nPhotons() > 0) nGathered = diffuseMap.gather(sp.P, gathered, nSearch, radius); 
			hp.radius2 = radius;
		}
	}

	state.userdata = o_udat;
	state.includeLights = oldIncludeLights;
	return;	
}


void SPPM::initializePPM(bool us_PM)
{
	const camera_t* camera = scene->getCamera();
	unsigned int resolution = camera->resX() * camera->resY();

	hitPoints.reserve(resolution);
	bound_t bBox = scene->getSceneBound(); // Now using Scene Bound, this could get a bigger initial radius, and need more tests

	// initialize SPPM statistics
	float initialRadius = ((bBox.longX() + bBox.longY() + bBox.longZ()) / 3.f) / ((camera->resX() + camera->resY()) / 2.0f) * 2.f ;
	initialRadius = std::min(initialRadius, 1.f); //Fix the overflow bug, but seems contain some bugs also, need more tests 
	for(int i = 0; i < resolution; i++)
	{
		HitPoint hp;

		hp.accPhotonFlux  = colorA_t(0.f);
		hp.accPhotonCount = 0;
		hp.radius2 = (initialRadius * initialFactor) * (initialRadius * initialFactor);

		hp.acclightsource = colorA_t(0.f);
		hp.constanthits = 0;
		hp.surfacehits = 0;

		hitPoints.push_back(hp);
	}
	
	if(bHashgrid) photonGrid.setParm(initialRadius*2.f, nPhotons, bBox);

}

integrator_t* SPPM::factory(paraMap_t &params, renderEnvironment_t &render)
{
	bool transpShad=false;
	bool pmIRE = false;
	int shadowDepth=5; //may used when integrate Direct Light
	int raydepth=5;
	int _passNum = 1000;
	int numPhotons = 500000;
	int bounces = 5;
	float times = 1.f;
	int searchNum = 100;

	params.getParam("transpShad", transpShad);
	params.getParam("shadowDepth", shadowDepth);
	params.getParam("raydepth", raydepth);
	params.getParam("photons", numPhotons);
	params.getParam("passNums", _passNum);
	params.getParam("bounces", bounces);
	params.getParam("times", times); // initial radius times

	params.getParam("searchNum", searchNum);
	params.getParam("pmIRE", pmIRE);

	SPPM* ite = new SPPM(numPhotons, _passNum,transpShad, shadowDepth);
	ite->rDepth = raydepth;
	ite->maxBounces = bounces;
	ite->initialFactor = times;

	ite->nSearch = searchNum; // under tests enable now
	ite->PM_IRE = pmIRE; 

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