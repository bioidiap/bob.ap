/**
 * @date Wed Jan 11:10:20 2013 +0200
 * @author Elie Khoury <Elie.Khoury@idiap.ch>
 * @author Laurent El Shafey <Laurent.El-Shafey@idiap.ch>
 * @author Pavel Korshunov <Pavel.Korshunov@idiap.ch>
 *
 * @brief Implement spectrogram
 *
 * Copyright (C) Idiap Research Institute, Martigny, Switzerland
 */


#ifndef BOB_AP_SPECTROGRAM_H
#define BOB_AP_SPECTROGRAM_H

#include <vector>
#include <stdexcept>
#include <blitz/array.h>
#include <boost/format.hpp>

#include <bob.sp/FFT1D.h>

#include "Energy.h"

namespace bob {
/**
 * \ingroup libap_api
 * @{
 *
 */
namespace ap {

/**
 * @brief This class implements an audio spectrogram extractor
 */
class Spectrogram: public Energy
{
  public:
    /**
     * @brief Constructor. Initializes working arrays
     */
    Spectrogram(const double sampling_frequency,
      const double win_length_ms=20., const double win_shift_ms=10.,
      const size_t n_filters=24, const double f_min=0.,
      const double f_max=8000., const double pre_emphasis_coeff=0.95,
      const bool mel_scale=true, const bool normalize_mean=true,
      const bool rect_filter=false, const bool inverse_filter=false,
      const bool normalize_spectrum=false, const bool ssfc_features=false,
      const bool scfc_features=false, const bool scmc_features=false);

    /**
     * @brief Copy Constructor
     */
    Spectrogram(const Spectrogram& other);

    /**
     * @brief Assignment operator
     */
    Spectrogram& operator=(const Spectrogram& other);

    /**
     * @brief Equal to
     */
    bool operator==(const Spectrogram& other) const;

    /**
     * @brief Not equal to
     */
    bool operator!=(const Spectrogram& other) const;

    /**
     * @brief Destructor
     */
    virtual ~Spectrogram();

    /**
     * @brief Gets the output shape for a given input/input length
     */
    virtual blitz::TinyVector<int,2> getShape(const size_t input_length) const;
    virtual blitz::TinyVector<int,2> getShape(const blitz::Array<double,1>& input) const;

    /**
     * @brief Computes the spectrogram
     */
    void operator()(const blitz::Array<double,1>& input, blitz::Array<double,2>& output);

    /**
     * @brief Returns the number of filters used in the filter bank.
     */
    size_t getNFilters() const
    { return m_n_filters; }
    /**
     * @brief Returns the frequency of the lowest triangular filter in the
     * filter bank
     */
    double getFMin() const
    { return m_f_min; }
    /**
     * @brief Returns the frequency of the highest triangular filter in the
     * filter bank
     */
    double getFMax() const
    { return m_f_max; }
    /**
     * @brief Tells whether the frequencies of the filters in the filter bank
     * are taken from the linear or the Mel scale
     */
    bool getMelScale() const
    { return m_mel_scale; }
    /**
     * @brief Tells whether the frequencies of the filters in the filter bank
     * are scaled using rectangular filter instead of Mel-scale or linear scale
     */
    bool getRectangularFilter() const
    { return m_rect_filter; }
    /**
     * @brief Tells whether to apply the filter in the inversed order, i.e., from high frequencies to low.
     */
    bool getInverseFilter() const
    {return m_inverse_filter;}
    /**
     * @brief Tells whether to normalize power spectrum.
     */
    bool getNormalizeSpectrum() const
    {return m_normalize_spectrum;}
    /**
     * @brief Tells whether SSFC features are being computed
     */
    bool getSSFCFeatures() const
    { return m_ssfc_features; }
    /**
     * @brief Tells whether SCFC features are being computed
     */
    bool getSCFCFeatures() const
    { return m_scfc_features; }
    /**
     * @brief Tells whether SCMC features are being computed
     */
    bool getSCMCFeatures() const
    { return m_scmc_features; }
    /**
     * @brief Returns the pre-emphasis coefficient.
     */
    double getPreEmphasisCoeff() const
    { return m_pre_emphasis_coeff; }
    /**
     * @brief Tells whether we used the energy or the square root of the energy
     */
    bool getEnergyFilter() const
    { return m_energy_filter; }
    /**
     * @brief Tells whether we used the log triangular filter or the triangular
     * filter
     */
    bool getLogFilter() const
    { return m_log_filter; }
    /**
     * @brief Tells whether we compute a spectrogram or energy bands
     */
    bool getEnergyBands() const
    { return m_energy_bands; }

    /**
     * @brief Sets the sampling frequency/frequency rate
     */
    virtual void setSamplingFrequency(const double sampling_frequency);
    /**
     * @brief Sets the window length in miliseconds
     */
    virtual void setWinLengthMs(const double win_length_ms);
    /**
     * @brief Sets the window shift in miliseconds
     */
    virtual void setWinShiftMs(const double win_shift_ms);

    /**
     * @brief Sets the number of filters used in the filter bank.
     */
    virtual void setNFilters(size_t n_filters);
    /**
     * @brief Sets the pre-emphasis coefficient. It should be a value in the
     * range [0,1].
     */
    virtual void setPreEmphasisCoeff(double pre_emphasis_coeff)
    {
      if (pre_emphasis_coeff < 0. || pre_emphasis_coeff > 1.) {
        boost::format m("the argument for `pre_emphasis_coeff' cannot take the value %f - the value must be in the interval [0.,1.]");
        m % pre_emphasis_coeff;
        throw std::runtime_error(m.str());
      }
      m_pre_emphasis_coeff = pre_emphasis_coeff;
    }
    /**
     * @brief Returns the frequency of the lowest triangular filter in the
     * filter bank
     */
    virtual void setFMin(double f_min);
    /**
     * @brief Returns the frequency of the highest triangular filter in the
     * filter bank
     */
    virtual void setFMax(double f_max);
    /**
     * @brief Sets whether the frequencies of the filters in the filter bank
     * are taken from the linear or the Mel scale
     */
    virtual void setMelScale(bool mel_scale);
    /**
     * @brief Sets whether the rectangular filter is used instead of the default linear
     */
    virtual void setRectangularFilter(bool rect_filter);
    /**
     * @brief Sets whether to apply the filter in the inversed order, i.e., from high frequencies to low.
     */
    virtual void setInverseFilter(bool inverse_filter);
    /**
     * @brief Sets whether to normalize the power spectrum of the signal.
     */
    virtual void setNormalizeSpectrum(bool normalize_spectrum)
    { m_normalize_spectrum = normalize_spectrum; }
    /**
     * @brief Set to true if you want to compute
     * Subband Spectral Flux Coefficients (SSFC), which measures
     * the frame-by-frame change in the power spectrum.
     */
    virtual void setSSFCFeatures(bool ssfc_features)
    { m_ssfc_features = ssfc_features; }
    /**
     * @brief Set to true if you want to compute
     * Spectral Centroid Frequency Coefficients (SCFC), which
     * capture detailed information about subbands similar to formant frequencies.
     */
    virtual void setSCFCFeatures(bool scfc_features);
    /**
     * @brief Set to true if you want to compute
     * Spectral Centroid Magnitude Coefficients (SCMC), which
     * capture detailed information about subbands similar to SCFC features.
     */
    virtual void setSCMCFeatures(bool scmc_features);

    /**
     * @brief Sets whether we used the energy or the square root of the energy
     */
    virtual void setEnergyFilter(bool energy_filter)
    { m_energy_filter = energy_filter; }
    /**
     * @brief Sets whether we used the log triangular filter or the triangular
     * filter
     */
    virtual void setLogFilter(bool log_filter)
    { m_log_filter = log_filter; }
    /**
     * @brief Sets whether we compute a spectrogram or energy bands
     */
    virtual void setEnergyBands(bool energy_bands)
    { m_energy_bands = energy_bands; }

  protected:
    /**
     * @brief Converts a frequency in Herz to the corresponding one in Mel
     */
    static double herzToMel(double f);
    /**
     * @brief Converts a frequency in Mel to the corresponding one in Herz
     */
    static double melToHerz(double f);
    /**
     * @brief Pre-emphasises the signal by applying the first order equation
     * \f$data_{n} := data_{n} − a*data_{n−1}\f$
     */
    void pre_emphasis(blitz::Array<double,1> &data, double& last_elem_prev_frame) const;
    /**
     * @brief Applies the Hamming window to the signal
     */
    void hammingWindow(blitz::Array<double,1> &data) const;

    /**
     * @brief Computes the power-spectrum of the FFT of the input frame
     */
    void powerSpectrumFFT(blitz::Array<double,1>& x);
    /**
     * @brief Applies the pre-computed filter bank (triangular or rectangular)
     */
    void filterBank(blitz::Array<double,1>& x);


    virtual void initWinLength();
    virtual void initWinSize();

    void initCacheHammingKernel();
    void initCacheFilterBank();

    /**
     * @brief Initialize the table m_p_index, which contains the indices of
     * the cut-off frequencies of the triangular filters.. It looks like:
     *
     *                      filter 2
     *                   <------------->
     *                filter 1           filter 4
     *             <----------->       <------------->
     *        | | | | | | | | | | | | | | | | | | | | | ..........
     *         0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9  ..........
     *             ^     ^     ^       ^             ^
     *             |     |     |       |             |
     *            t[0]   |    t[2]     |           t[4]
     *                  t[1]          t[3]
     *
     */
    void initCachePIndex();
    void initCacheFilters();

    size_t m_n_filters;
    double m_f_min;
    double m_f_max;
    double m_pre_emphasis_coeff;
    bool m_mel_scale;
    bool m_rect_filter;
    bool m_inverse_filter;
    bool m_normalize_spectrum;
    bool m_ssfc_features; //this flag should be true for SSFC features computation
    bool m_scfc_features; //this flag should be true for SCFC features computation
    bool m_scmc_features; //this flag should be true for SCMC features computation
    double m_fb_out_floor;
    bool m_energy_filter;
    bool m_log_filter;
    bool m_energy_bands;
    double m_log_fb_out_floor;

    blitz::Array<double,1> m_hamming_kernel;
    blitz::Array<double,1> m_p_index;
    std::vector<blitz::Array<double,1> > m_filter_bank;
    std::vector<blitz::Array<double,1> > m_filter_weights;
    bob::sp::FFT1D m_fft;

    mutable blitz::Array<std::complex<double>,1> m_cache_frame_c1;
    mutable blitz::Array<std::complex<double>,1> m_cache_frame_c2;
    mutable blitz::Array<double,1> m_cache_filters;
};

}
}

#endif /* BOB_AP_SPECTROGRAM_H */
