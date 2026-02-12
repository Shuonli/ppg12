#ifndef MBDPILEUPHELPER_H
#define MBDPILEUPHELPER_H

#include <cmath>
#include <algorithm>

/**
 * @brief Pileup detection based on MBD timing information
 * 
 * Adapted from PileupRejector.cc - uses timing spread across MBD PMTs
 * to identify pileup events.
 */

enum class PileupCutStrength {
    COMFORT,    // Uses product of RMS times (less strict)
    STRICT,     // Uses average of RMS times
    DRACONIAN   // Uses average of RMS times with tighter cut
};

struct MbdPileupResult {
    float chargesum = 0;
    float prodsigma = 0;   // product of north/south RMS times
    float avgsigma = 0;    // average of north/south RMS times
    float maxsigma = 0;    // max of north/south RMS times
    float proddelta = 0;   // product of (max-min) time differences
    float avgdelta = 0;    // average of (max-min) time differences
    float maxdelta = 0;    // max of (max-min) time differences
    float npmt[2] = {0, 0};  // number of PMTs passing cuts per side
    bool valid = false;    // true if enough PMTs on at least one side
};

/**
 * @brief Calculate pileup metrics from MBD PMT timing arrays
 * 
 * @param northt   Array of north MBD PMT times (64 channels)
 * @param southt   Array of south MBD PMT times (64 channels)
 * @param northq   Array of north MBD PMT charges (64 channels)
 * @param southq   Array of south MBD PMT charges (64 channels)
 * @param hitcut   Minimum number of hits per side to compute timing (default: 2)
 * @param time_cut Maximum |time| to accept PMT (default: 25 ns)
 * @param charge_cut Minimum charge to accept PMT (default: 0.4)
 * @return MbdPileupResult containing all computed metrics
 */
template<typename ArrayType>
MbdPileupResult calculateMbdPileupMetrics(
    const ArrayType& northt,
    const ArrayType& southt,
    const ArrayType& northq,
    const ArrayType& southq,
    int hitcut = 2,
    float time_cut = 25.0f,
    float charge_cut = 0.4f)
{
    MbdPileupResult result;
    
    float maxtime[2] = {-100, -100};
    float mintime[2] = {100, 100};
    float rmstime[2] = {0, 0};
    float sumtime[2] = {0, 0};
    
    // Process north side (side 0)
    const int nchannels = 64;
    for (int i = 0; i < nchannels; i++)
    {
        float charge = northq[i];
        float time = northt[i];
        
        if (std::fabs(time) < time_cut && charge > charge_cut)
        {
            result.chargesum += charge;
            result.npmt[0]++;
            if (time > maxtime[0]) maxtime[0] = time;
            if (time < mintime[0]) mintime[0] = time;
            rmstime[0] += time * time;
            sumtime[0] += time;
        }
    }
    
    // Process south side (side 1)
    for (int i = 0; i < nchannels; i++)
    {
        float charge = southq[i];
        float time = southt[i];
        
        if (std::fabs(time) < time_cut && charge > charge_cut)
        {
            result.chargesum += charge;
            result.npmt[1]++;
            if (time > maxtime[1]) maxtime[1] = time;
            if (time < mintime[1]) mintime[1] = time;
            rmstime[1] += time * time;
            sumtime[1] += time;
        }
    }
    
    // Check if we have enough hits
    if (result.npmt[0] < 1 && result.npmt[1] < 1)
    {
        result.valid = false;
        return result;
    }
    result.valid = true;
    
    // Calculate RMS for north side
    if (result.npmt[0] >= hitcut)
    {
        rmstime[0] /= result.npmt[0];
        sumtime[0] /= result.npmt[0];
        float variance = rmstime[0] - sumtime[0] * sumtime[0];
        rmstime[0] = variance > 0 ? std::sqrt(variance) : 0;
    }
    else
    {
        maxtime[0] = 0;
        mintime[0] = 0;
        rmstime[0] = 0;
        sumtime[0] = 0;
    }
    
    // Calculate RMS for south side
    if (result.npmt[1] >= hitcut)
    {
        rmstime[1] /= result.npmt[1];
        sumtime[1] /= result.npmt[1];
        float variance = rmstime[1] - sumtime[1] * sumtime[1];
        rmstime[1] = variance > 0 ? std::sqrt(variance) : 0;
    }
    else
    {
        maxtime[1] = 0;
        mintime[1] = 0;
        rmstime[1] = 0;
        sumtime[1] = 0;
    }
    
    // Calculate final metrics
    result.prodsigma = rmstime[0] * rmstime[1];
    result.avgsigma = (rmstime[0] + rmstime[1]) / 2.0f;
    result.maxsigma = std::max(rmstime[0], rmstime[1]);
    
    result.proddelta = (maxtime[1] - mintime[1]) * (maxtime[0] - mintime[0]);
    result.avgdelta = ((maxtime[1] - mintime[1]) + (maxtime[0] - mintime[0])) / 2.0f;
    result.maxdelta = std::max(maxtime[1] - mintime[1], maxtime[0] - mintime[0]);
    
    return result;
}

/**
 * @brief Determine if event is pileup based on calculated metrics
 * 
 * @param result   The pileup metrics from calculateMbdPileupMetrics
 * @param strength The cut strength to apply
 * @param comfort_cut  Cut value for COMFORT mode (default: 1.0)
 * @param strict_cut   Cut value for STRICT mode (default: 2.0)
 * @param draconian_cut Cut value for DRACONIAN mode (default: 1.5)
 * @return true if event is identified as pileup
 */
inline bool isMbdPileup(
    const MbdPileupResult& result,
    PileupCutStrength strength = PileupCutStrength::STRICT,
    float comfort_cut = 1.0f,
    float strict_cut = 2.0f,
    float draconian_cut = 1.5f)
{
    if (!result.valid)
    {
        return false;
    }
    
    switch (strength)
    {
        case PileupCutStrength::COMFORT:
            return result.prodsigma >= comfort_cut;
        case PileupCutStrength::STRICT:
            return result.avgsigma >= strict_cut;
        case PileupCutStrength::DRACONIAN:
            return result.avgsigma >= draconian_cut;
        default:
            return false;
    }
}

/**
 * @brief Convenience function to check pileup in one call
 * 
 * @param northt   Array of north MBD PMT times (64 channels)
 * @param southt   Array of south MBD PMT times (64 channels)
 * @param northq   Array of north MBD PMT charges (64 channels)
 * @param southq   Array of south MBD PMT charges (64 channels)
 * @param strength The cut strength to apply
 * @param hitcut   Minimum number of hits per side
 * @return true if event is identified as pileup
 */
template<typename ArrayType>
bool checkMbdPileup(
    const ArrayType& northt,
    const ArrayType& southt,
    const ArrayType& northq,
    const ArrayType& southq,
    PileupCutStrength strength = PileupCutStrength::STRICT,
    int hitcut = 2)
{
    MbdPileupResult result = calculateMbdPileupMetrics(northt, southt, northq, southq, hitcut);
    return isMbdPileup(result, strength);
}

#endif // MBDPILEUPHELPER_H
