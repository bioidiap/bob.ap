/**
 * @date Wed Jan 11:09:30 2013 +0200
 * @author Elie Khoury <Elie.Khoury@idiap.ch>
 * @author Laurent El Shafey <Laurent.El-Shafey@idiap.ch>
 * @author Pavel Korshunov <Pavel.Korshunov@idiap.ch>
 *
 *
 * Copyright (C) Idiap Research Institute, Martigny, Switzerland
 */

#include <bob.ap/Spectrogram.h>
#include <bob.core/assert.h>
#include <bob.core/cast.h>

bob::ap::Spectrogram::Spectrogram(const double sampling_frequency,
    const double win_length_ms, const double win_shift_ms,
    const size_t n_filters, const double f_min, const double f_max,
    const double pre_emphasis_coeff, const bool mel_scale,
    const bool normalize_mean, const bool rect_filter,
    const bool inverse_filter, const bool normalize_spectrum,
    const bool ssfc_features, const bool scfc_features, const bool scmc_features):
  bob::ap::Energy(sampling_frequency, win_length_ms, win_shift_ms, normalize_mean),
  m_n_filters(n_filters), m_f_min(f_min), m_f_max(f_max),
  m_pre_emphasis_coeff(pre_emphasis_coeff), m_mel_scale(mel_scale),
  m_rect_filter(rect_filter), m_inverse_filter(inverse_filter),
  m_normalize_spectrum(normalize_spectrum), m_ssfc_features(ssfc_features),
  m_scfc_features(scfc_features), m_scmc_features(scmc_features),
  m_fb_out_floor(1.), m_energy_filter(false), m_log_filter(true),
  m_energy_bands(false), m_fft()
{
  // Check pre-emphasis coefficient
  if (pre_emphasis_coeff < 0. || pre_emphasis_coeff > 1.) {
    boost::format m("the argument for `pre_emphasis_coeff' cannot take the value %f - the value must be in the interval [0.,1.]");
    m % pre_emphasis_coeff;
    throw std::runtime_error(m.str());
  }

  // Initialization
  initWinLength();
  initWinShift();

  // Initializes logarithm of flooring values
  // pavel - allow computing log for a very small value
  m_fb_out_floor = std::numeric_limits<double>::epsilon();
  m_log_fb_out_floor = log(m_fb_out_floor);

  m_cache_filters.resize(m_n_filters);
}

bob::ap::Spectrogram::Spectrogram(const Spectrogram& other):
  bob::ap::Energy(other), m_n_filters(other.m_n_filters),
  m_f_min(other.m_f_min), m_f_max(other.m_f_max),
  m_pre_emphasis_coeff(other.m_pre_emphasis_coeff),
  m_mel_scale(other.m_mel_scale),
  m_rect_filter(other.m_rect_filter),
  m_inverse_filter(other.m_inverse_filter),
  m_normalize_spectrum(other.m_normalize_spectrum),
  m_ssfc_features(other.m_ssfc_features),
  m_scfc_features(other.m_scfc_features),
  m_scmc_features(other.m_scmc_features),
  m_fb_out_floor(other.m_fb_out_floor),
  m_energy_filter(other.m_energy_filter),
  m_log_filter(other.m_log_filter),
  m_energy_bands(other.m_energy_bands),
  m_fft(other.m_fft)
{
  // Initialization
  initWinLength();
  initWinShift();

  // Initializes logarithm of flooring values
  m_log_fb_out_floor = log(m_fb_out_floor);

  m_cache_filters.resize(m_n_filters);
}

bob::ap::Spectrogram& bob::ap::Spectrogram::operator=(const bob::ap::Spectrogram& other)
{
  if (this != &other)
  {
    bob::ap::Energy::operator=(other);
    m_n_filters = other.m_n_filters;
    m_f_min = other.m_f_min;
    m_f_max = other.m_f_max;
    m_pre_emphasis_coeff = other.m_pre_emphasis_coeff;
    m_mel_scale = other.m_mel_scale;
    m_rect_filter = other.m_rect_filter;
    m_inverse_filter = other.m_inverse_filter;
    m_normalize_spectrum = other.m_normalize_spectrum;
    m_ssfc_features = other.m_ssfc_features;
    m_scfc_features = other.m_scfc_features;
    m_scmc_features = other.m_scmc_features;
    m_fb_out_floor = other.m_fb_out_floor;
    m_energy_filter = other.m_energy_filter;
    m_log_filter = other.m_log_filter;
    m_energy_bands = other.m_energy_bands;
    m_fft = other.m_fft;

    // Initialization
    initWinLength();
    initWinShift();

    // Initializes logarithm of flooring values
    m_log_fb_out_floor = log(m_fb_out_floor);

    m_cache_filters.resize(m_n_filters);
  }
  return *this;
}

bool bob::ap::Spectrogram::operator==(const bob::ap::Spectrogram& other) const
{
  return (bob::ap::Energy::operator==(other) &&
          m_n_filters == other.m_n_filters && m_f_min == other.m_f_min &&
          m_f_max == other.m_f_max &&
          m_pre_emphasis_coeff == other.m_pre_emphasis_coeff &&
          m_mel_scale == other.m_mel_scale &&
          m_rect_filter == other.m_rect_filter &&
          m_normalize_spectrum == other.m_normalize_spectrum &&
          m_inverse_filter == other.m_inverse_filter &&
          m_ssfc_features == other.m_ssfc_features &&
          m_scfc_features == other.m_scfc_features &&
          m_scmc_features == other.m_scmc_features &&
          m_fb_out_floor == other.m_fb_out_floor &&
          m_energy_filter == other.m_energy_filter &&
          m_log_filter == other.m_log_filter &&
          m_energy_bands == other.m_energy_bands);
}

bool bob::ap::Spectrogram::operator!=(const bob::ap::Spectrogram& other) const
{
  return !(this->operator==(other));
}

bob::ap::Spectrogram::~Spectrogram()
{
}

blitz::TinyVector<int,2>
bob::ap::Spectrogram::getShape(const size_t input_size) const
{
  // Res will contain the number of frames x the dimension of the feature vector
  blitz::TinyVector<int,2> res;

  // 1. Number of frames
  res(0) = 1+((input_size-m_win_length)/m_win_shift);

  // 2. Dimension of the feature vector
  res(1) = (m_energy_bands? m_n_filters : m_win_size/2 + 1);

  return res;
}

blitz::TinyVector<int,2>
bob::ap::Spectrogram::getShape(const blitz::Array<double,1>& input) const
{
  return getShape(input.extent(0));
}

void bob::ap::Spectrogram::setSamplingFrequency(const double sampling_frequency)
{
  bob::ap::Energy::setSamplingFrequency(sampling_frequency);
  initWinLength();
  initWinShift();
}

void bob::ap::Spectrogram::setWinLengthMs(const double win_length_ms)
{
  bob::ap::Energy::setWinLengthMs(win_length_ms);
  initWinLength();
}

void bob::ap::Spectrogram::setWinShiftMs(const double win_shift_ms)
{
  bob::ap::Energy::setWinShiftMs(win_shift_ms);
  initWinShift();
}

void bob::ap::Spectrogram::setNFilters(size_t n_filters)
{
  m_n_filters = n_filters;
  m_cache_filters.resize(m_n_filters);
  initCacheFilterBank();
}

void bob::ap::Spectrogram::setFMin(double f_min)
{
  m_f_min = f_min;
  initCacheFilterBank();
}

void bob::ap::Spectrogram::setFMax(double f_max)
{
  m_f_max = f_max;
  initCacheFilterBank();
}

void bob::ap::Spectrogram::setMelScale(bool mel_scale)
{
  m_mel_scale = mel_scale;
  initCacheFilterBank();
}

void bob::ap::Spectrogram::setInverseFilter(bool inverse_filter)
{
  m_inverse_filter = inverse_filter;
  initCacheFilterBank();
}

void bob::ap::Spectrogram::setRectangularFilter(bool rect_filter)
{
  m_rect_filter = rect_filter;
  initCacheFilterBank();
}

void bob::ap::Spectrogram::setSCFCFeatures(bool scfc_features)
{
  m_scfc_features = scfc_features;
  initCacheFilters();
}

void bob::ap::Spectrogram::setSCMCFeatures(bool scmc_features)
{
  m_scmc_features = scmc_features;
  initCacheFilters();
}

double bob::ap::Spectrogram::herzToMel(double f)
{
  return (2595.*log10(1+f/700.));
}

double bob::ap::Spectrogram::melToHerz(double f)
{
  return ((double)(700.*(pow(10,f/2595.)-1)));
}

void bob::ap::Spectrogram::initCacheHammingKernel()
{
  // Hamming Window initialization
  m_hamming_kernel.resize(m_win_length);
  double cst = 2*M_PI/(double)(m_win_length-1);
  blitz::firstIndex i;
  m_hamming_kernel = 0.54-0.46*blitz::cos(i*cst);
}

void bob::ap::Spectrogram::initCacheFilterBank()
{
  initCachePIndex();
  initCacheFilters();
}

void bob::ap::Spectrogram::initCachePIndex()
{
  // Computes the indices for the triangular filter bank
  m_p_index.resize(m_n_filters+2);
  // 'Mel' frequency decomposition (for MFCC)
  if (m_mel_scale)
  {
    double m_max = herzToMel(m_f_max);
    double m_min = herzToMel(m_f_min);
    for (int i=0; i<(int)m_n_filters+2; ++i) {
      double alpha = i/ (double)(m_n_filters+1);
      double f = melToHerz(m_min * (1-alpha) + m_max * alpha);
//      if (m_inverse_filter)
//        m_p_index(m_n_filters+1-i)=(int)round(m_win_size * (m_f_max - f) / m_sampling_frequency); // from max to min
//      else {
        double factor = f / m_sampling_frequency;
      // pavel - use double values instead of integers for better precision
      // m_p_index(i)=(int)round(m_win_size * factor); // normal Mel-filter from min to max
        m_p_index(i)=m_win_size * factor; // normal Mel-filter from min to max
//      }
    }
  }
  else
  // Linear frequency decomposition (for LFCC or RFCC)
  {
    const double cst_a = (m_win_size/m_sampling_frequency) * (m_f_max-m_f_min)/(double)(m_n_filters+1);
    const double cst_b = (m_win_size/m_sampling_frequency) * m_f_min;
    for (int i=0; i<(int)m_n_filters+2; ++i) {
      // pavel - temporarily change round() to floor() for compatibility
//      m_p_index(i) = (int)floor(cst_a * i + cst_b);
      // use double values instead of integers for better precision
      m_p_index(i) = cst_a * i + cst_b;
    }
  }
}

void bob::ap::Spectrogram::initCacheFilters()
{
  // Creates the Triangular filter bank
  m_filter_bank.clear();
  m_filter_weights.clear();
  blitz::firstIndex ii;
  for (int i=0; i<(int)m_n_filters; ++i)
  {
    // Integer indices of the boundary of the triangular filter in the
    // Fourier domain
    //make sure the left border of the interval is not-included, except for the first i=0
    int li = (int)floor(m_p_index(i)+1);
    int mi = (int)floor(m_p_index(i+1));
    int ri = (int)floor(m_p_index(i+2));
    if (i == 0 || (ri-li == 0)) //first elem or left and last elements are the same
      li--;
    blitz::Array<double,1> filt(ri-li+1);

    // if we have a rectangular filter, we do not need to care about slices
    if (m_rect_filter) {
      // Rectangular filter with value=1. for all indexes
      filt = 1.;
    }
    else {

      // Fill in the left slice of the triangular filter
      blitz::Array<double,1> filt_p1(filt(blitz::Range(0,mi-li)));
      double denom = m_p_index(i+1)-m_p_index(i);
      filt_p1 = (ii+li-m_p_index(i)) / denom;
      // Fill in the right slice of the triangular filter
      blitz::Array<double,1> filt_p2(filt(blitz::Range(mi-li+1,ri-li)));
      denom = m_p_index(i+2)-m_p_index(i+1);
      filt_p2 = (m_p_index(i+2)-(ii+mi+1)) / denom;
    }
    // Append filter into the filterbank vector
    m_filter_bank.push_back(filt);
    // init the weights used for some features during filtering
    // create an array of size equal to the current range
    blitz::Array<double,1> weights(ri-li+1);
    weights = 1.;
    // compute normalized frequencies for SSFC or SCMC features computation
    if (m_scfc_features || m_scmc_features) {
      // li=m_p_index(i) is the first adjusted frequency of the current range
      // as a normalizer, we use the maximum possible adjusted frequency value
      // so, we divide by (2*m_win_size/m_sampling_frequency) * (m_f_max-m_f_min)
      // in a case when m_f_min=0 and m_f_max is 8000 (maximum), the denominator will become 512 - the max possible range
      weights = 1.0*(li + ii) / ((2*m_win_size/m_sampling_frequency) * (m_f_max-m_f_min)); //make sure it's double value
    }
    // Append current weights to the weights vector
    m_filter_weights.push_back(weights);
  }
}

void bob::ap::Spectrogram::initWinLength()
{
  bob::ap::Energy::initWinLength();
  initCacheHammingKernel();
  initCacheFilterBank();
}

void bob::ap::Spectrogram::initWinSize()
{
  bob::ap::Energy::initWinSize();
  m_fft.setLength(m_win_size);
  m_cache_frame_c1.resize(m_win_size);
  m_cache_frame_c2.resize(m_win_size);
}

void bob::ap::Spectrogram::pre_emphasis(blitz::Array<double,1> &data, double& last_elem_prev_frame) const
{
  if (m_pre_emphasis_coeff != 0.)
  {
    // remember the last element of the frame that will be needed
    // to compute emphasis correctly in the next frame
    // the element is at the end of the m_win_shift, since next frame starts at m_win_shift position
    double last_element = data((int)m_win_shift-1);

    // Pre-emphasise the signal by applying the first order equation
    // \f$data_{n} := data_{n} − a*data_{n−1}\f$
    blitz::Range r0((int)m_win_length-2,0,-1);
    blitz::Range r1((int)m_win_length-1,1,-1);
    data(r1) -= m_pre_emphasis_coeff * data(r0); // Apply first order equation
    // pavel - remove first element update for consistency with Matlab
//    data(0) *= 1. - m_pre_emphasis_coeff; // Update first element
    // pavel - use the the last remembered element of the previous frame to update the fist element of this frame
    data(0) -= m_pre_emphasis_coeff * last_elem_prev_frame; // Update first elementm_last_element
    // remember the last element of this frame for the next frame
    last_elem_prev_frame = last_element;
  }
}

void bob::ap::Spectrogram::hammingWindow(blitz::Array<double,1> &data) const
{
  blitz::Range r(0,(int)m_win_length-1);
  data(r) *= m_hamming_kernel;
}

void bob::ap::Spectrogram::powerSpectrumFFT(blitz::Array<double,1>& x)
{
  // this function, effectivelly, modifies only the first half of the input array 'x'
  // first, copy the input into a complex container
  m_cache_frame_c1 = bob::core::array::cast<std::complex<double> >(x);
  // Apply the FFT and store the results in complex m_cache_frame_c2
  m_fft(m_cache_frame_c1, m_cache_frame_c2);

  // Take the the power spectrum of the first part of the output of the FFT
  // This range (half of the original array) is what we work on from now on
  blitz::Range r(0,(int)m_win_size/2);
  // point (do not copy) to the first half of the original input
  blitz::Array<double,1> x_half(x(r));
  // point to the first half of the result of FFT
  blitz::Array<std::complex<double>,1> complex_half(m_cache_frame_c2(r));
  // basically, abs(of complex array) gives us the magnitude of the power spectrum
  x_half = blitz::abs(complex_half);
  if (m_energy_filter) // Energy is basically magnitude in power of 2
    x_half = blitz::pow2(x_half);
  if (m_normalize_spectrum) {// Normalize power spectrum, if we need the normalized value
    //    x_half -= blitz::mean(x_half); // this is an older implementation but it may give similar results
    // it should be maximum of the FFT over the whole signal but we have access only to one frame at a time
    double max_fft = blitz::max(x_half);
    if (max_fft > std::numeric_limits<double>::epsilon())
      x_half /= max_fft;
  }
}

void bob::ap::Spectrogram::filterBank(blitz::Array<double,1>& x)
{
  // first, create an array that we operate on (half of the power spectrum)
  blitz::Array<double,1> x_half(x(blitz::Range(0,(int)m_win_size/2)));

  // apply a reversed filters to the data (e.g., IMFCC-inversed of MFCC). The easiest way is to
  // reverse the order of the data, while keeping everything, including filter banks, intact
  if (m_inverse_filter)
    x_half.reverseSelf(blitz::firstDim); //reverse only the first half of the power spectrum, since the other half is ignored anyway

  for (int i=0; i<(int)m_n_filters; ++i)
  {

    // pavel - ensure that we use each frequiency only once!
    // it means we use interval with non-inclusive left border
    int first_fr = (int)floor(m_p_index(i)+1);
    int last_fr = (int)floor(m_p_index(i+2));
    // except for the very first interval, when we start with first (or zeros) m_p_index index.
    // or when first and last index are the same
    if (i == 0 || first_fr == last_fr)
      first_fr--;
    // take the pre-computed range of frequencies corresponding to the bank 'i'
    blitz::Range slice_range = blitz::Range(first_fr, last_fr);
    // take the slice of data corresponding to those frequencies
    blitz::Array<double,1> data_slice(x_half(slice_range));

    // apply pre-computed bank filter on the data (which should be power spectrum) and
    // multiply by pre-computed m_filter_weights[] (equal to 1 except for SCFC or SCMC features)
    // carefull not to modify data_slice, since it is just a range of the input x vector
    // and we have overlapping windows here, so we should not change the input x
    double res = blitz::sum(m_filter_weights[i] * data_slice * m_filter_bank[i]);

    // for SSFC features, divide the result by the sum of filtered magnitude of the power spectrum
    if (m_scfc_features)
      res /= blitz::sum(data_slice * m_filter_bank[i]);
    // for SCMC features, divide the result by the sum of normalized frequencies
    if (m_scmc_features)
      res /= blitz::sum(m_filter_weights[i]);

    // store the result in m_cache_filters, which will be later used in the computation pipeline
    if (m_log_filter)
      // Take the log of the result, i.e., compute log triangular filter bank
      m_cache_filters(i)= (res < m_fb_out_floor ? m_log_fb_out_floor : log(res));
    else
      m_cache_filters(i) = res;
  }

}

void bob::ap::Spectrogram::operator()(const blitz::Array<double,1>& input,
  blitz::Array<double,2>& spectrogram_matrix)
{
  // Get expected dimensionality of output array
  blitz::TinyVector<int,2> spectrogram_shape = bob::ap::Spectrogram::getShape(input);
  // Check dimensionality of output array
  bob::core::array::assertSameShape(spectrogram_matrix, spectrogram_shape);
  int n_frames=spectrogram_shape(0);

  // Computes the center of the cut-off frequencies
  blitz::Range r1 = blitz::Range(0,m_win_size/2);
  if (m_energy_bands)
    r1 = blitz::Range(0,m_n_filters-1);

  double last_frame_elem=0;

  for (int i=0; i<n_frames; ++i)
  {
    // Extract and normalize frame
    extractNormalizeFrame(input, i, m_cache_frame_d);

    // Apply pre-emphasis
    pre_emphasis(m_cache_frame_d, last_frame_elem);

    // Apply the Hamming window
    hammingWindow(m_cache_frame_d);

    // Take the power spectrum of the first part of the FFT
    powerSpectrumFFT(m_cache_frame_d);

    // Filter with the triangular filter bank (either in linear or Mel domain)
    if (m_energy_bands)
      filterBank(m_cache_frame_d);

    blitz::Array<double,1> spec_matrix_row(spectrogram_matrix(i,r1));
    if (m_energy_bands)
      spec_matrix_row = m_cache_filters(r1);
    else
      spec_matrix_row = m_cache_frame_d(r1);
  }
}


