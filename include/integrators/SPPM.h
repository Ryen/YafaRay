#ifndef Y_SPPM_H
#define Y_SPPM_H

#include <yafray_config.h>
#include <core_api/environment.h>
#include <core_api/material.h>
#include <core_api/background.h>
#include <core_api/light.h>
#include <core_api/imagefilm.h>
#include <core_api/camera.h>
#include <core_api/mcintegrator.h>
#include <yafraycore/photon.h>
#include <yafraycore/monitor.h>
#include <yafraycore/ccthreads.h>
#include <yafraycore/timer.h>
#include <yafraycore/spectrum.h>
#include <utilities/sample_utils.h>
#include <utilities/mcqmc.h>
#include <yafraycore/scr_halton.h>

__BEGIN_YAFRAY

typedef struct shared_statistics // per-pixel amount
{
	float radius2; // square radius,
	unsigned long long accPhotonCount;
	colorA_t accPhotonFlux; // accumulated flux

	//for constant_color usage
	colorA_t acclightsource; // accumulated radiance from light source & emiting material
	unsigned int surfacehits; //not used now
	unsigned int constanthits; // same as above
}HitPoint;

typedef struct GInfo
{
	unsigned long photonCount;
	colorA_t photonFlux;
	colorA_t constantRandiance;

	GInfo(){photonCount = 0; photonFlux = colorA_t(0.f); constantRandiance = colorA_t(0.f);}

	GInfo & operator +=(const GInfo &g)
	{ 
		photonCount += g.photonCount; 
		photonFlux += g.photonFlux;
		constantRandiance += g.constantRandiance;
		return (*this);
	}
}GatherInfo;


class YAFRAYPLUGIN_EXPORT SPPM: public mcIntegrator_t
{
	public:
	SPPM(unsigned int dPhotons, int passnum);
	~SPPM();
	virtual bool render(imageFilm_t *imageFilm);
	/*! render a pass; only required by the default implementation of render() */
	virtual bool renderPass(int samples, int offset, bool adaptive);
	/*! render a tile; only required by default implementation of render() */
	virtual bool renderTile(renderArea_t &a, int n_samples, int offset, bool adaptive, int threadID);
	virtual bool preprocess(); //not used for now
	virtual void prePass(int samples, int offset, bool adaptive); // do a  photon pass
	virtual colorA_t integrate(renderState_t &state, diffRay_t &ray/*, sampler_t &sam*/) const; // not used now
	static integrator_t* factory(paraMap_t &params, renderEnvironment_t &render);

	void initializePPM(bool us_PM);
	GatherInfo traceGatherRay(renderState_t &state, diffRay_t &ray, HitPoint &hp); //based on integrate method to do the gatering trace, need double-check deadly.
	void traceIRERay(renderState_t &state, diffRay_t &ray, HitPoint &hp); //based on integrate method

	protected:
	photonMap_t diffuseMap,causticMap; // photonmap
	std::string settings;
	pdf1D_t *lightPowerD;
	unsigned int nPhotons; //photon number to scatter
	int nSearch;// now used to do inital radius estimate
	int passNum; // the progressive pass number
	float initialFactor; // used to time the initial radius
	unsigned int totalnPhotons; // amount used to normalize photon energy
	//unsigned int totalnCausPhotons;
	bool PM_IRE; //use PM for initial radius estimate

	Halton hal2,hal3,hal5,hal7, hal8,hal9,hal10;
	unsigned int nRefined; // Debug info: Refined pixel per pass

	std::vector<shared_statistics>hitPoints; // per-pixel refine data

};

__END_YAFRAY

#endif // SPPM