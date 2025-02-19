// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

// PhysLight code contributed by Anders Langlands and Luca Fascione
// Copyright (c) 2020, Weta Digital, Ltd.
// SPDX-License-Identifier: Apache-2.0

// P3D updates

#include <mitsuba/render/estimator.h>
NAMESPACE_BEGIN(mitsuba)

MTS_VARIANT std::unique_ptr<Estimator<Float,Spectrum>> Estimator<Float,Spectrum>::Create(const std::string &name) {

    std::unique_ptr<Estimator<Float,Spectrum>> estimator;

    // TODO: Later use of paramset maybe..
    if (name == "mean")
        estimator = std::make_unique<MeanEstimator<Float,Spectrum>>(name);
    else if (name == "mon")
        estimator = std::make_unique<MONEstimator<Float,Spectrum>>(name);
    // else if (name == "gini-mon")
    //     estimator = std::make_unique<GiniMONEstimator>(name);
    // else if (name == "gini-binary-mon")
    //     estimator = std::make_unique<GiniBinaryMONEstimator>(name);
    // else if (name == "gini-partial-mon")
    //     estimator = std::make_unique<GiniPartialMONEstimator>(name);
    // else if (name == "gini-dmon") // <---
    //     estimator = std::make_unique<GiniDistMONEstimator>(name);
    // else if (name == "gini-partial-dmon")
    //     estimator = std::make_unique<GiniDistPartialMONEstimator>(name);
    // else if (name == "pakmon")
    //     estimator = std::make_unique<PakMONEstimator>(name);
    // else if (name == "mean-or-mon")
    //     estimator = std::make_unique<MeanOrMONEstimator>(name);
    // else if (name == "abmm")
    //     estimator = std::make_unique<ABMMEstimator>(name);
    // else if (name == "gabmm")
    //     estimator = std::make_unique<GABMMEstimator>(name);
    else {
        printf("%s: estimator type unknown. Use of default: mean", name.c_str());
        estimator = std::make_unique<MeanEstimator<Float,Spectrum>>(name);
    }

    if (!estimator)
        printf("%s: unable to create estimator.", name.c_str());

    return estimator;
}

MTS_VARIANT std::string Estimator<Float,Spectrum>::ToString() const {
    return name + "Estimator";
}

MTS_VARIANT void MeanEstimator<Float,Spectrum>::Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const
{
    
    weightSum = 0.;

    // get each weightSum of pixelMoN
    for (int j = 0; j < pixelWindow.windowSize; j++) {
        weightSum += pixelWindow.buffers[j].weightSum;
    }

    // based on channel numbers
    for (int i = 0; i < 3; i++) {

        // loop over pixels (used as means storage) for computing real channel value
        rgb[i] = 0.;
        splatRGB[i] = 0.;

        for (int j = 0; j < pixelWindow.windowSize; j++) {
            rgb[i] += pixelWindow.buffers[j].rgbSum[i];
            splatRGB[i] = splatRGB[i] + pixelWindow.buffers[j].splatRGB[i];
        }
    }
};


MTS_VARIANT void MONEstimator<Float,Spectrum>::Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const
{
    
    weightSum = 0.;

    // based on channel numbers
    for (int i = 0; i < 3; i++) {

        // store channel information
        std::vector<Float> cvalues;
        std::vector<Float> weightsSum;
        std::vector<double> csplats;

        for (int j = 0; j < pixelWindow.windowSize; j++) {
            cvalues.push_back(pixelWindow.buffers[j].rgbSum[i]);
            // per channel management (but weight can be different depending of median buffer)
            weightsSum.push_back(pixelWindow.buffers[j].weightSum);
            csplats.push_back(pixelWindow.buffers[j].splatRGB[i]);
        }

        // temp storage in order to sort values
        std::vector<Float> means(cvalues);
        std::vector<int> sortedIndices = means.sort();

        Float monWeight, monMean = 0.;
        double monSplat = 0;

        // compute median from means
        // find associated weightsum index and use it
        // Classical MON
        if (nbuffers % 2 == 1){
            unsigned unsortedIndex = sortedIndices[int(nbuffers/2)];

            monMean = cvalues[unsortedIndex];
            monWeight = weightsSum[unsortedIndex];
            monSplat = csplats[unsortedIndex];
        }
        else{
            int k_mean = int(nbuffers/2);
            unsigned firstIndex = sortedIndices[k_mean - 1];
            unsigned secondIndex = sortedIndices[k_mean];

            monMean = (cvalues[firstIndex] + cvalues[secondIndex]) / 2;
            monWeight = (weightsSum[firstIndex] + weightsSum[secondIndex]) / 2;
            monSplat = (csplats[firstIndex]  + csplats[secondIndex]) / 2;
        }

        // store channel information
        weightSum += monWeight;
        rgb[i] = monMean;
        splatRGB[i] = monSplat;
    }

    // divide per number of channel the weightSum
    weightSum /= 3;
};


// MTS_VARIANT void ABMMEstimator::Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const {
//     this->Estimate(pixelWindow, rgb, weightSum, splatRGB, 1.);
// };

// MTS_VARIANT void ABMMEstimator::Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB, Float alpha) const
// {
    
//     weightSum = 0.;

//     // get each weightSum of pixelMoN
//     for (int j = 0; j < pixelWindow.windowSize; j++) {
//         weightSum += pixelWindow.buffers[j].weightSum;
//     }

//     // based on channel numbers
//     for (int i = 0; i < 3; i++) {

//         // loop over pixels (used as means storage) for computing real channel value
//         rgb[i] = 0.;
//         splatRGB[i] = 0.;
        
//         Float sum = 0.;
//         Float squaredSum = 0;
//         Float cubicSum = 0;

//         for (int j = 0; j < pixelWindow.windowSize; j++) {
//             rgb[i] += pixelWindow.buffers[j].rgbSum[i];
//             sum += pixelWindow.buffers[j].rgbSum[i];
//             splatRGB[i] = splatRGB[i] + pixelWindow.buffers[j].splatRGB[i];
//             squaredSum += pixelWindow.buffers[j].squaredSum[i];
//             cubicSum += pixelWindow.buffers[j].cubicSum[i];
//         }

//         if (sum > 0.) {
//             // computation of different moments
//             int n = pixelWindow.nsamples; // using nsamples

//             Float mean = sum / n; // M1
//             Float onlineM2 = (squaredSum / n) - (mean * mean); // M2
//             Float stdTheta = std::sqrt(onlineM2);

//             Float onlineM3 = (cubicSum - 3 * mean * squaredSum) / n + 2 * (mean * mean * mean);
//             Float onlineSkew = onlineM3 / std::pow(onlineM2, 1.5);

//             // std::cout << "------------------------" << std::endl;
//             // std::cout << "Weight is: " << n << std::endl;
//             // std::cout << "Sum is: " << sum << std::endl;
//             // std::cout << "SquaredSum is: " << squaredSum << std::endl;
//             // std::cout << "CubicSum is: " << cubicSum << std::endl;
//             // std::cout << "Mean is: " << mean << std::endl;
//             // std::cout << "M2 is: " << onlineM2 << std::endl;
//             // std::cout << "M3 is: " << onlineM3 << std::endl;
//             // std::cout << "Skew is: " << onlineSkew << std::endl;
//             Float aBMM = mean - ((1. / 3.) * (stdTheta / (n * alpha + 2))) * onlineSkew;
//             // std::cout << "Estimation is: " << aBBM << std::endl;
//             rgb[i] = aBMM * n;
//         }
//     }
// };

// void GABMMEstimator::Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const
// {
//     // for each channel number
//     Float giniSum = 0;

//     for (int i = 0; i < 3; i++) {

//         std::vector<Float> cvalues;

//         for (int j = 0; j < pixelWindow.windowSize; j++) {
//             cvalues.push_back(pixelWindow.buffers[j].rgbSum[i] / pixelWindow.buffers[j].weightSum);
//         }

//         giniSum += this->getGini(cvalues);
//     }

//     Float giniMean = giniSum / 3.;

//     // use of gini ean for estimate
//     aBMMEstimator->Estimate(pixelWindow, rgb, weightSum, splatRGB, 1. - giniMean);
// };

// void PakMONEstimator::Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const
// {
//     weightSum = 0.;

//     // based on channel numbers
//     for (int i = 0; i < 3; i++) {

//         // store channel information
//         std::vector<Float> cvalues;
//         std::vector<Float> weightsSum;
//         std::vector<double> csplats;

//         for (int j = 0; j < pixelWindow.windowSize; j++) {
//             cvalues.push_back(pixelWindow.buffers[j].rgbSum[i]);
//             // per channel management (but weight can be different depending of median buffer)
//             weightsSum.push_back(pixelWindow.buffers[j].weightSum);
//             csplats.push_back(pixelWindow.buffers[j].splatRGB[i]);
//         }

//         // temp storage in order to sort values
//         std::vector<Float> means(cvalues);
//         std::vector<int> sortedIndices = means.sort();

//         // sum storage
//         Float meansSum = 0;

//         // PakMON expected output
//         Float weight, mean = 0.;
//         double csplat = 0;

//         // by default classical MON values
//         if (nbuffers % 2 == 1){
//             unsigned unsortedIndex = sortedIndices[int(nbuffers/2)];

//             mean = cvalues[unsortedIndex];
//             weight = weightsSum[unsortedIndex];
//             csplat = csplats[unsortedIndex];
//         }
//         else{
//             int k_mean = int(nbuffers/2);
//             unsigned firstIndex = sortedIndices[k_mean - 1];
//             unsigned secondIndex = sortedIndices[k_mean];

//             mean = (cvalues[firstIndex] + cvalues[secondIndex]) / 2;
//             weight = (weightsSum[firstIndex] + weightsSum[secondIndex]) / 2;
//             csplat = (csplats[firstIndex]  + csplats[secondIndex]) / 2;
//         }

//         for (int j = 0; j < cvalues.size(); j++)
//             meansSum += cvalues[j];

//         Float currentMean = meansSum / cvalues.size();
            
//         // compute variance distance evolution
//         std::vector<Float> distances;

//         for (int j = 2; j < means.size(); j++) {
            
//             // use of sorted means in order to compute variance evolution by step 1
//             // compute variance of each elements
//             Float var = 0;
            
//             // use of previously sorted means
//             for(int k = 0; k < j; k++)
//             {
//                 var += (means[k] - currentMean) * (means[k] - currentMean);
//             }
//             var /= (j + 1);

//             // add new 
//             distances.push_back(var);
//         }

//         // use of variance evolution and compute entropy
//         Float distancesEntropy = getEntropy(distances);

//         // Computation of PakMON using \alpha and \rho value
//         unsigned middleIndex = int(nbuffers / 2);

//         // alpha and rho automatically set value
//         Float alpha = distancesEntropy;

//         if (alpha < 0.000000001) {
//             alpha = 0.000000001;
//         }
//         int rho = (int)(middleIndex * distancesEntropy) - (int)(nbuffers * 0.15); // try using avoid 30% (total) of current nbuffers

//         unsigned lowerIndex = 0;
//         unsigned higherIndex = 0;
    
//         // get current lower and higher index 
//         if (nbuffers % 2 == 0) {
            
//             lowerIndex = middleIndex - 1;
//             higherIndex = middleIndex;
            
//         } else {
//             lowerIndex = middleIndex;
//             higherIndex = middleIndex;
//         }

//         // use of sorted means and relative sorted indices
//         for (int j = 1; j < rho + 1; j++) {

//             // current neighbor multiple factor
//             Float multFactor = pow(alpha, j);

//             // add left and right neighbor contribution
//             mean += means[lowerIndex - j] * multFactor;
//             mean += means[higherIndex + j] * multFactor;
            
//             // weighting contribution to take in account
//             // use of this index to retrieve the associated weightsSum
//             weight += weightsSum[sortedIndices[lowerIndex]] * multFactor;
//             weight += weightsSum[sortedIndices[higherIndex]] * multFactor;

//             csplat += csplats[sortedIndices[lowerIndex]] * multFactor;
//             csplat += csplats[sortedIndices[higherIndex]] * multFactor;
//         }

//         // store channel information
//         weightSum += weight;
//         rgb[i] = mean;
//         splatRGB[i] = csplat;
//     }

//     // divide per number of channel the weightSum
//     weightSum /= 3;
// };

// Float PakMONEstimator::getEntropy(std::vector<Float> values) const {

//     // computation of squared values
//     Float sumEigenValues = 0;
//     std::vector<Float> eigenValues(values.size());

//     for (int i = 0; i < values.size(); i++) {
//         Float sum = values[i] * values[i];
//         eigenValues[i] = sum;
//         sumEigenValues += sum;
//     }

//     // normalization the squared values
//     std::vector<Float> v(values.size());

//     for (int i = 0; i < values.size(); i++) {
//         // add of epsilon value
//         v[i] = eigenValues[i] / (sumEigenValues + 0.000000000001);
//     }

//     // computation of entropy
//     Float entropy = 0;

//     for (int i = 0; i < values.size(); i++) {
//         if (v[i] > 0) {
//             entropy += v[i] * log(v[i]);
//         }
//     }

//     entropy *= -1;

//     entropy /= log(values.size());

//     return entropy;
// };


// void AlphaMONEstimator::Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const
// {
//    this->Estimate(pixelWindow, rgb, weightSum, splatRGB, 0.5); // default use of 0.5 confidence
// };

// void AlphaMONEstimator::Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB, Float alpha) const
// {
//     weightSum = 0.;

//     // based on channel numbers
//     for (int i = 0; i < 3; i++) {

//         // store channel information
//         std::vector<Float> cvalues;
//         std::vector<Float> weightsSum;
//         std::vector<double> csplats;

//         for (int j = 0; j < pixelWindow.windowSize; j++) {
//             cvalues.push_back(pixelWindow.buffers[j].rgbSum[i]);
//             // per channel management (but weight can be different depending of median buffer)
//             weightsSum.push_back(pixelWindow.buffers[j].weightSum);
//             csplats.push_back(pixelWindow.buffers[j].splatRGB[i]);
//         }

//         // temp storage in order to sort values
//         std::vector<Float> means(cvalues);
//         std::vector<int> sortedIndices = means.sort();

//         // sum storage
//         Float meansSum = 0;

//         // PakMON expected output
//         Float weight, mean = 0.;
//         double csplat = 0;

//         // by default classical MON values
//         if (nbuffers % 2 == 1){
//             unsigned unsortedIndex = sortedIndices[int(nbuffers/2)];

//             mean = cvalues[unsortedIndex];
//             weight = weightsSum[unsortedIndex];
//             csplat = csplats[unsortedIndex];
//         }
//         else{
//             int k_mean = int(nbuffers/2);
//             unsigned firstIndex = sortedIndices[k_mean - 1];
//             unsigned secondIndex = sortedIndices[k_mean];

//             mean = (cvalues[firstIndex] + cvalues[secondIndex]) / 2;
//             weight = (weightsSum[firstIndex] + weightsSum[secondIndex]) / 2;
//             csplat = (csplats[firstIndex]  + csplats[secondIndex]) / 2;
//         }

//         // Computation of PakMON using \alpha and \rho value
//         unsigned middleIndex = int(nbuffers / 2);

//         unsigned lowerIndex = 0;
//         unsigned higherIndex = 0;
//         unsigned until = middleIndex;
    
//         // get current lower and higher index 
//         if (nbuffers % 2 == 0) {
            
//             lowerIndex = middleIndex - 1;
//             higherIndex = middleIndex;
//             until = middleIndex - 1;
            
//         } else {
//             lowerIndex = middleIndex;
//             higherIndex = middleIndex;
//         }

//         // use of sorted means and relative sorted indices
//         for (int j = 1; j < until + 1; j++) {

//             // current neighbor multiple factor
//             Float multFactor = pow(alpha, j);

//             // add left and right neighbor contribution
//             mean += means[lowerIndex - j] * multFactor;
//             mean += means[higherIndex + j] * multFactor;
            
//             // weighting contribution to take in account
//             // use of this index to retrieve the associated weightsSum
//             weight += weightsSum[sortedIndices[lowerIndex]] * multFactor;
//             weight += weightsSum[sortedIndices[higherIndex]] * multFactor;

//             csplat += csplats[sortedIndices[lowerIndex]] * multFactor;
//             csplat += csplats[sortedIndices[higherIndex]] * multFactor;
//         }

//         // store channel information
//         weightSum += weight;
//         rgb[i] = mean;
//         splatRGB[i] = csplat;
//     }

//     // divide per number of channel the weightSum
//     weightSum /= 3;
// };

// void AlphaDistMONEstimator::Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const
// {
//    this->Estimate(pixelWindow, rgb, weightSum, splatRGB, 0.5); // default use of 0.5 confidence
// };

// void AlphaDistMONEstimator::Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB, Float alpha) const
// {
//     weightSum = 0.;

//     // based on channel numbers
//     for (int i = 0; i < 3; i++) {

//         // store channel information
//         std::vector<Float> cvalues;
//         std::vector<Float> weightsSum;
//         std::vector<double> csplats;

//         for (int j = 0; j < pixelWindow.windowSize; j++) {
//             cvalues.push_back(pixelWindow.buffers[j].rgbSum[i]);
//             // per channel management (but weight can be different depending of median buffer)
//             weightsSum.push_back(pixelWindow.buffers[j].weightSum);
//             csplats.push_back(pixelWindow.buffers[j].splatRGB[i]);
//         }

//         // temp storage in order to sort values
//         std::vector<Float> means(cvalues);
//         std::vector<int> sortedIndices = means.sort();

//         // sum storage
//         Float meansSum = 0;

//         // PakMON expected output
//         Float weight, mean = 0.;
//         double csplat = 0;

//         // by default classical MON values
//         if (nbuffers % 2 == 1){
//             unsigned unsortedIndex = sortedIndices[int(nbuffers/2)];

//             mean = cvalues[unsortedIndex];
//             weight = weightsSum[unsortedIndex];
//             csplat = csplats[unsortedIndex];
//         }
//         else{
//             int k_mean = int(nbuffers/2);
//             unsigned firstIndex = sortedIndices[k_mean - 1];
//             unsigned secondIndex = sortedIndices[k_mean];

//             mean = (cvalues[firstIndex] + cvalues[secondIndex]) / 2;
//             weight = (weightsSum[firstIndex] + weightsSum[secondIndex]) / 2;
//             csplat = (csplats[firstIndex]  + csplats[secondIndex]) / 2;
//         }

//         // Computation of PakMON using \alpha and \rho value
//         unsigned middleIndex = int(nbuffers / 2);

//         unsigned lowerIndex = 0;
//         unsigned higherIndex = 0;
//         unsigned until = int(middleIndex * alpha);
    
//         // get current lower and higher index 
//         if (nbuffers % 2 == 0) {
            
//             lowerIndex = middleIndex - 1;
//             higherIndex = middleIndex;
//             until = middleIndex - 1;
            
//         } else {
//             lowerIndex = middleIndex;
//             higherIndex = middleIndex;
//         }

//         // use of sorted means and relative sorted indices
//         for (int j = 1; j < until + 1; j++) {

//             // add left and right neighbor contribution
//             mean += means[lowerIndex - j];
//             mean += means[higherIndex + j];
            
//             // weighting contribution to take in account
//             // use of this index to retrieve the associated weightsSum
//             weight += weightsSum[sortedIndices[lowerIndex]];
//             weight += weightsSum[sortedIndices[higherIndex]];

//             csplat += csplats[sortedIndices[lowerIndex]];
//             csplat += csplats[sortedIndices[higherIndex]];
//         }

//         // store channel information
//         weightSum += weight;
//         rgb[i] = mean;
//         splatRGB[i] = csplat;
//     }

//     // divide per number of channel the weightSum
//     weightSum /= 3;
// };

// void MeanOrMONEstimator::Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const
// {
//     // Check use of mon or use of mean based on IC
//     Float meanICSum = 0.;
//     Float monICSum = 0.;

//     // based on channel numbers
//     for (int i = 0; i < 3; i++) {

//         // loop over pixels (used as means storage) for computing real channel value
//         std::vector<Float> cvalues;
//         std::vector<Float> csquared;
//         std::vector<Float> weightsSum;
//         std::vector<double> csplats;

//         for (int j = 0; j < pixelWindow.windowSize; j++) {
//             cvalues.push_back(pixelWindow.buffers[j].rgbSum[i]);
//             csquared.push_back(pixelWindow.buffers[j].squaredSum[i]);
//             weightsSum.push_back(pixelWindow.buffers[j].weightSum);
//             csplats.push_back(pixelWindow.buffers[j].splatRGB[i]);
//         }

//         std::vector<Float> means(cvalues); // copy of current channel values

//         // depending of confidence interval of predictor, compute the required value
//         // need to compute mean value, hence all values are needed
//         Float meanSquaredValues = 0.;
        
//         Float meanWeight, meanMean = 0.;
//         double meanSplat = 0;

//         for (int j = 0; j < nbuffers; j++) {
            
//             meanMean += means[j];
//             meanSquaredValues += csquared[j];
//             meanWeight += weightsSum[j];
//             meanSplat += csplats[j];
//         }

//         /////////////////////////////////
//         // confidence interval of Mean //
//         /////////////////////////////////
//         Float currentMeanMean = meanMean / meanWeight;

//         // use of weight as number of samples
//         Float meanStdValue = (meanSquaredValues / meanWeight) - (currentMeanMean * currentMeanMean);
//         Float meanIC = (1.96 * meanStdValue) / std::sqrt(meanWeight);

//         /////////////////////////////////
//         // confidence interval of MON  //
//         /////////////////////////////////

//         // sort means vector and get sorted indices as output
//         Float monSquaredValue = 0;

//         Float monWeight, monMean = 0.;
//         double monSplat = 0;

//         std::vector<int> sortedIndices = means.sort();

//         // compute median from means
//         // find associated weightsum index and use it
//         // Classical MON
//         if (nbuffers % 2 == 1){
//             unsigned unsortedIndex = sortedIndices[int(nbuffers/2)];

//             monWeight = weightsSum[unsortedIndex];
//             monMean = cvalues[unsortedIndex];
//             monSquaredValue = csquared[unsortedIndex];
//             monSplat = csplats[unsortedIndex];
//         }
//         else{
//             int k_mean = int(nbuffers/2);
//             unsigned firstIndex = sortedIndices[k_mean - 1];
//             unsigned secondIndex = sortedIndices[k_mean];

//             monWeight = (weightsSum[firstIndex] + weightsSum[secondIndex]) / 2;
//             monMean = (cvalues[firstIndex] + cvalues[secondIndex]) / 2;
//             monSquaredValue = (csquared[firstIndex] + csquared[secondIndex]) / 2;
//             monSplat = (csplats[firstIndex]  + csplats[secondIndex]) / 2;
//         }

//         Float currentMonMean = monMean / monWeight;
//         Float monStdValue = (monSquaredValue / monWeight) - (currentMonMean * currentMonMean);
//         Float monIC = (1.96 * monStdValue) / std::sqrt(monWeight);

//         meanICSum += meanIC;
//         monICSum += monIC;
//     }

//     // check use of mean or mon
//     if ((meanICSum / 3)  <= (monICSum / 3)) {
//         meanEstimator->Estimate(pixelWindow, rgb, weightSum, splatRGB);
//     } else {
//         monEstimator->Estimate(pixelWindow, rgb, weightSum, splatRGB);
//     }
// };


// Float GiniMONEstimator::getEntropy(std::vector<Float> values) const {

//     // computation of squared values
//     Float sumEigenValues = 0;
//     std::vector<Float> eigenValues(values.size());

//     for (int i = 0; i < values.size(); i++) {
//         Float sum = values[i] * values[i];
//         eigenValues[i] = sum;
//         sumEigenValues += sum;
//     }

//     // normalization the squared values
//     std::vector<Float> v(values.size());

//     for (int i = 0; i < values.size(); i++) {
//         // add of epsilon value
//         v[i] = eigenValues[i] / (sumEigenValues + 0.000000000001);
//     }

//     // computation of entropy
//     Float entropy = 0;

//     for (int i = 0; i < values.size(); i++) {
//         if (v[i] > 0) {
//             entropy += v[i] * log(v[i]);
//         }
//     }

//     entropy *= -1;

//     entropy /= log(values.size());

//     return entropy;
// };

// Float GiniMONEstimator::getGini(std::vector<Float> values) const {

//     // get indices of array
//     int n = values.size();
//     Float arraySum = 0;
//     Float indexArraySum = 0;

//     Float minValue = Infinity;

//     // get min value
//     for (int i = 0; i < n; i++)
//         if (values[i] < minValue)
//             minValue = values[i];

//     // need to sort obtained values
//     std::sort(values.begin(), values.end());

//     // avoid 0 value and store index
//     for (int i = 0; i < n; i++) {

//         // avoid negative value
//         if (minValue < 0)
//             values[i] -= minValue; 

//         values[i] += 0.00000000001; // epsilon value
//         arraySum += values[i];
//         indexArraySum += (2 * (i + 1) - n - 1) * values[i];
//     }

//     return indexArraySum / (n * arraySum);
// }

// void GiniMONEstimator::Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const
// {
//     // for each channel number
//     Float giniSum = 0;

//     for (int i = 0; i < 3; i++) {

//         std::vector<Float> cvalues;

//         for (int j = 0; j < pixelWindow.windowSize; j++) {
//             cvalues.push_back(pixelWindow.buffers[j].rgbSum[i] / pixelWindow.buffers[j].weightSum);
//         }

//         giniSum += this->getGini(cvalues);
//     }

//     Float giniMean = giniSum / 3.;

//     alphaMoNEstimator->Estimate(pixelWindow, rgb, weightSum, splatRGB, 1 - giniMean);
// };

// void GiniDistMONEstimator::Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const
// {
//     // for each channel number
//     Float giniSum = 0;

//     for (int i = 0; i < 3; i++) {

//         std::vector<Float> cvalues;

//         for (int j = 0; j < pixelWindow.windowSize; j++) {
//             cvalues.push_back(pixelWindow.buffers[j].rgbSum[i] / pixelWindow.buffers[j].weightSum);
//         }

//         giniSum += this->getGini(cvalues);
//     }

//     Float giniMean = giniSum / 3.;

//     alphaDistMoNEstimator->Estimate(pixelWindow, rgb, weightSum, splatRGB, 1 - giniMean);
// };

// Float GiniDistMONEstimator::getGini(std::vector<Float> values) const {

//     // get indices of array
//     int n = values.size();
//     Float arraySum = 0;
//     Float indexArraySum = 0;

//     Float minValue = Infinity;

//     // get min value
//     for (int i = 0; i < n; i++)
//         if (values[i] < minValue)
//             minValue = values[i];

//     // need to sort obtained values
//     std::sort(values.begin(), values.end());

//     // avoid 0 value and store index
//     for (int i = 0; i < n; i++) {

//         // avoid negative value
//         if (minValue < 0)
//             values[i] -= minValue; 

//         values[i] += 0.00000000001; // epsilon value
//         arraySum += values[i];
//         indexArraySum += (2 * (i + 1) - n - 1) * values[i];
//     }

//     return indexArraySum / (n * arraySum);
// }

// void GiniBinaryMONEstimator::Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const
// {
//     // for each channel number
//     Float giniSum = 0;

//     for (int i = 0; i < 3; i++) {

//         std::vector<Float> cvalues;

//         for (int j = 0; j < pixelWindow.windowSize; j++) {
//             cvalues.push_back(pixelWindow.buffers[j].rgbSum[i] / pixelWindow.buffers[j].weightSum);
//         }

//         giniSum += this->getGini(cvalues);
//     }

//     Float giniMean = giniSum / 3.;

//     if (giniMean < 0.25)
//         alphaMoNEstimator->Estimate(pixelWindow, rgb, weightSum, splatRGB, 1.);
//     else
//         alphaMoNEstimator->Estimate(pixelWindow, rgb, weightSum, splatRGB, 0.);
// };


// void GiniPartialMONEstimator::Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const
// {
//     // for each channel number
//     Float giniSum = 0;

//     for (int i = 0; i < 3; i++) {

//         std::vector<Float> cvalues;

//         for (int j = 0; j < pixelWindow.windowSize; j++) {
//             cvalues.push_back(pixelWindow.buffers[j].rgbSum[i] / pixelWindow.buffers[j].weightSum);
//         }

//         giniSum += this->getGini(cvalues);
//     }

//     Float giniMean = giniSum / 3.;

//     if (giniMean < 0.25)
//         alphaMoNEstimator->Estimate(pixelWindow, rgb, weightSum, splatRGB, 1.);
//     else
//         alphaMoNEstimator->Estimate(pixelWindow, rgb, weightSum, splatRGB, 0.5);
// };

// void GiniDistPartialMONEstimator::Estimate(const PixelWindow &pixelWindow,  Spectrum &rgb, Float &weightSum, double* splatRGB) const
// {
//     // for each channel number
//     Float giniSum = 0;

//     for (int i = 0; i < 3; i++) {

//         std::vector<Float> cvalues;

//         for (int j = 0; j < pixelWindow.windowSize; j++) {
//             cvalues.push_back(pixelWindow.buffers[j].rgbSum[i] / pixelWindow.buffers[j].weightSum);
//         }

//         giniSum += this->getGini(cvalues);
//     }

//     Float giniMean = giniSum / 3.;

//     if (giniMean < 0.25)
//         alphaDistMoNEstimator->Estimate(pixelWindow, rgb, weightSum, splatRGB, 1.);
//     else
//         alphaDistMoNEstimator->Estimate(pixelWindow, rgb, weightSum, splatRGB, 0.5);
// };

NAMESPACE_END(mitsuba)
