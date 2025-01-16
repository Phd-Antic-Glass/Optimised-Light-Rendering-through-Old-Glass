// Code from p3d_PBRT-v4
// Jerome Buisine
// ported to Mitsuba2 by Quentin Huan

#pragma once
#include <mitsuba/mitsuba.h>
#include <mitsuba/core/logger.h>
#include <mitsuba/core/object.h>
#include <mitsuba/core/vector.h>

#include <mitsuba/core/rfilter.h>
#include <mitsuba/core/vector.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/fwd.h>


#include <atomic>
#include <string>
#include <vector>
#include <math.h>
NAMESPACE_BEGIN(mitsuba)

// P3D update  mon parameter (if value of 1, no kmon use, hence classical mean)
const int nbuffers = 11;

// Need to externalise PixelWindow declaration for RGBFilm
// P3D Updates
struct PixelBuffer {
    PixelBuffer() = default;

    double rgbSum[3] = {0., 0., 0.};
    double squaredSum[3] = {0., 0., 0.};
    double cubicSum[3] = {0., 0., 0.};
    double splatRGB[3];
    double weightSum = 0.;

    void Clear() {

        for (int i = 0; i < 3; i++) {
            rgbSum[i] = 0.;
            squaredSum[i] = 0.;
            cubicSum[i] = 0.;
            splatRGB[i] = 0.;
        }

        weightSum = 0.;
    }
};

// P3D Updates
struct PixelWindow {
    PixelWindow() = default;

    PixelBuffer buffers[nbuffers];

    int windowSize = nbuffers; // number of buffers clusters
    int index = 0; // keep track of index used
    int nsamples = 0; // keep track of nsamples;
    bool filled = false;
};

// Base Estimator class MTS_EXPORT_RENDER
template <typename Float, typename Spectrum>
class MTS_EXPORT_RENDER MTS_EXPORT_RENDER Estimator : public Object {

public:
    MTS_IMPORT_TYPES()

        ~Estimator() {};
        
        static std::unique_ptr<Estimator> Create(const std::string &name); 

         
        virtual void Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const = 0;

        std::string ToString() const;
    MTS_DECLARE_CLASS()

    protected:

        Estimator(const std::string &name) : name(name) {};
       
        std::string name;
};

// // approximated Bayesian Median of Means Estimator class MTS_EXPORT_RENDER MTS_EXPORT_RENDER
// template <typename Float, typename Spectrum>
// class MTS_EXPORT_RENDER ABMMEstimator : public Estimator<Float,Spectrum> {

// public:
//     MTS_IMPORT_TYPES()

//         ABMMEstimator(const std::string &name) : Estimator(name) {}; 

         
//         void Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const;

         
//         void Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB, Float alpha) const;
//     MTS_DECLARE_CLASS()

// };

// Mean Estimator class MTS_EXPORT_RENDER
template <typename Float, typename Spectrum>
class MTS_EXPORT_RENDER MeanEstimator : public Estimator<Float,Spectrum> {
public:
    MTS_IMPORT_BASE(Estimator)
    MTS_IMPORT_TYPES()

        MeanEstimator(const std::string &name) : Estimator<Float,Spectrum>(name) {}; 

         
        void Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const;
    MTS_DECLARE_CLASS()

};

// MON Estimator class MTS_EXPORT_RENDER
// Median of meaNs: use of median value from available mean buffers
template <typename Float, typename Spectrum>
class MTS_EXPORT_RENDER MONEstimator : public Estimator<Float,Spectrum> {

public:
    MTS_IMPORT_BASE(Estimator)
    MTS_IMPORT_TYPES()

        MONEstimator(const std::string &name) : Estimator<Float,Spectrum>(name) {}; 

         
        void Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const;
    MTS_DECLARE_CLASS()

};

// // AlphaMON Estimator class MTS_EXPORT_RENDER
// // Median of meaNs: use of median value from available mean buffers
// // use of an alpha criterion for convergence
// template <typename Float, typename Spectrum>
// class MTS_EXPORT_RENDER AlphaMONEstimator : public Estimator {

// public:
//     MTS_IMPORT_TYPES()

//         AlphaMONEstimator(const std::string &name) : Estimator(name) {}; 

         
//         void Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const;

         
//         void Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB, Float alpha) const;
//     MTS_DECLARE_CLASS()

// };

// // AlphaDistMON Estimator class MTS_EXPORT_RENDER
// // Median of meaNs: use of median value from available mean buffers
// // use of an alpha criterion for convergence using whole package
// template <typename Float, typename Spectrum>
// class MTS_EXPORT_RENDER AlphaDistMONEstimator : public Estimator {

// public:
//     MTS_IMPORT_TYPES()

//         AlphaDistMONEstimator(const std::string &name) : Estimator(name) {}; 

         
//         void Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const;

         
//         void Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB, Float alpha) const;
//     MTS_DECLARE_CLASS()

// };

// // GiniMON Estimator class MTS_EXPORT_RENDER
// // Median of meaNs: use of median value from available mean buffers
// // Use of Gini in order to well use \alpha criterion
// template <typename Float, typename Spectrum>
// class MTS_EXPORT_RENDER GiniDistMONEstimator : public Estimator {

// public:
//     MTS_IMPORT_TYPES()

//         GiniDistMONEstimator(const std::string &name) : Estimator(name) {

//             // default alpha value
//             alphaDistMoNEstimator = std::make_unique<AlphaDistMONEstimator>("admon");
//         }; 

         
//         void Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const;
//     MTS_DECLARE_CLASS()

//     protected:
//         std::unique_ptr<AlphaDistMONEstimator> alphaDistMoNEstimator;

         
//         Float getGini(pstd::vector<Float> values) const;
// };

// template <typename Float, typename Spectrum>
// class MTS_EXPORT_RENDER GiniDistPartialMONEstimator : public GiniDistMONEstimator {

//     MTS_IMPORT_TYPES()
// public:

//         GiniDistPartialMONEstimator(const std::string &name) : GiniDistMONEstimator(name) {};

         
//         void Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const;
// };

// // GiniMON Estimator class MTS_EXPORT_RENDER
// // Median of meaNs: use of median value from available mean buffers
// // Use of Gini in order to well use \alpha criterion
// template <typename Float, typename Spectrum>
// class MTS_EXPORT_RENDER GiniMONEstimator : public Estimator {

// public:
//     MTS_IMPORT_TYPES()

//         GiniMONEstimator(const std::string &name) : Estimator(name) {

//             // default alpha value
//             alphaMoNEstimator = std::make_unique<AlphaMONEstimator>("amon");
//         }; 

         
//         void Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const;
//     MTS_DECLARE_CLASS()

//     protected:
//         std::unique_ptr<AlphaMONEstimator> alphaMoNEstimator;

         
//         Float getEntropy(pstd::vector<Float> values) const;

         
//         Float getGini(pstd::vector<Float> values) const;
// };

// // GiniaBMM Estimator class MTS_EXPORT_RENDER
// // Median of meaNs: use of median value from available mean buffers
// // Use of Gini in order to well use \alpha criterion
// template <typename Float, typename Spectrum>
// class MTS_EXPORT_RENDER GABMMEstimator : public GiniMONEstimator {

// public:
//     MTS_IMPORT_TYPES()

//         GABMMEstimator(const std::string &name) : GiniMONEstimator(name) {

//             // default alpha value
//             aBMMEstimator = std::make_unique<ABMMEstimator>("abmm");
//         }; 

         
//         void Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const;
//     MTS_DECLARE_CLASS()

//     protected:
//         std::unique_ptr<ABMMEstimator> aBMMEstimator;
// };

// template <typename Float, typename Spectrum>
// class MTS_EXPORT_RENDER GiniBinaryMONEstimator : public GiniMONEstimator {

// public:
//     MTS_IMPORT_TYPES()

//         GiniBinaryMONEstimator(const std::string &name) : GiniMONEstimator(name) {};

         
//         void Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const;
//     MTS_DECLARE_CLASS()

// };

// template <typename Float, typename Spectrum>
// class MTS_EXPORT_RENDER GiniPartialMONEstimator : public GiniMONEstimator {

// public:
//     MTS_IMPORT_TYPES()

//         GiniPartialMONEstimator(const std::string &name) : GiniMONEstimator(name) {};

         
//         void Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const;
//     MTS_DECLARE_CLASS()

// };


// // PakMON Estimation class MTS_EXPORT_RENDER
// // Based of MON Estimator but with confidence criteria of median's neighborhood buffers
// template <typename Float, typename Spectrum>
// class MTS_EXPORT_RENDER PakMONEstimator : public Estimator {

// public:
//     MTS_IMPORT_TYPES()

//         PakMONEstimator(const std::string &name) : Estimator(name) {}; 

         
//         void Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const;
//     MTS_DECLARE_CLASS()

//     private:
         
//         Float getEntropy(pstd::vector<Float> values) const;

// };

// // Mean or MON estimator
// // Based on Confidence Interval criteria of pixel
// // Final estimator is chosen (mean or MON)
// template <typename Float, typename Spectrum>
// class MTS_EXPORT_RENDER MeanOrMONEstimator : public Estimator {

// public:
//     MTS_IMPORT_TYPES()

//         MeanOrMONEstimator(const std::string &name) : Estimator(name) {
//             meanEstimator = std::make_unique<MeanEstimator>("mean");
//             monEstimator = std::make_unique<MONEstimator>("mon");
//         }; 

         
//         void Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const;
//     MTS_DECLARE_CLASS()

//     private:
//         std::unique_ptr<MeanEstimator> meanEstimator;
//         std::unique_ptr<MONEstimator> monEstimator;
// };

NAMESPACE_END(mitsuba)
