/**
 * @date Wed Jan 11:09:30 2013 +0200
 * @author Elie Khoury <Elie.Khoury@idiap.ch>
 * @author Laurent El Shafey <Laurent.El-Shafey@idiap.ch>
 * @author Pavel Korshunov <Pavel.Korshunov@idiap.ch>
*
 * Copyright (C) Idiap Research Institute, Martigny, Switzerland
 */

#include <stdexcept>

#include <bob.ap/FrameExtractor.h>
#include <bob.core/check.h>

bob::ap::FrameExtractor::FrameExtractor(const double sampling_frequency,
    const double win_length_ms, const double win_shift_ms,
    const bool normalize_mean):
  m_sampling_frequency(sampling_frequency), m_win_length_ms(win_length_ms),
  m_win_shift_ms(win_shift_ms), m_normalize_mean(normalize_mean)
{
  // Initialization
  initWinLength();
  initWinShift();
  initMaxRange();
}

bob::ap::FrameExtractor::FrameExtractor(const FrameExtractor& other):
  m_sampling_frequency(other.m_sampling_frequency),
  m_win_length_ms(other.m_win_length_ms),
  m_win_shift_ms(other.m_win_shift_ms),
  m_normalize_mean(other.m_normalize_mean)
{
  // Initialization
  initWinLength();
  initWinShift();
  initMaxRange();
}

bob::ap::FrameExtractor::~FrameExtractor()
{
}

bob::ap::FrameExtractor& bob::ap::FrameExtractor::operator=(const bob::ap::FrameExtractor& other)
{
  if (this != &other)
  {
    m_sampling_frequency = other.m_sampling_frequency;
    m_win_length_ms = other.m_win_length_ms;
    m_win_shift_ms = other.m_win_shift_ms;
    m_normalize_mean = other.m_normalize_mean;

    // Initialization
    initWinLength();
    initWinShift();
    initMaxRange();
  }
  return *this;
}

bool bob::ap::FrameExtractor::operator==(const bob::ap::FrameExtractor& other) const
{
  return (m_sampling_frequency == other.m_sampling_frequency &&
      m_win_length_ms == other.m_win_length_ms &&
      m_win_shift_ms == other.m_win_shift_ms &&
      m_normalize_mean == other.m_normalize_mean);
}

bool bob::ap::FrameExtractor::operator!=(const bob::ap::FrameExtractor& other) const
{
  return !(this->operator==(other));
}

void bob::ap::FrameExtractor::setSamplingFrequency(const double sampling_frequency)
{
  m_sampling_frequency = sampling_frequency;
  initWinLength();
  initWinShift();
  initMaxRange();
}

void bob::ap::FrameExtractor::setWinLengthMs(const double win_length_ms)
{
  m_win_length_ms = win_length_ms;
  initWinLength();
}

void bob::ap::FrameExtractor::setWinShiftMs(const double win_shift_ms)
{
  m_win_shift_ms = win_shift_ms;
  initWinShift();
}

void bob::ap::FrameExtractor::initWinLength()
{
  m_win_length = (size_t)(m_sampling_frequency * m_win_length_ms / 1000);
  if (m_win_length == 0)
    throw std::runtime_error("The length of the window is 0. You should use a larger sampling rate or window length in miliseconds");
  initWinSize();

}

void bob::ap::FrameExtractor::initWinShift()
{
  m_win_shift = (size_t)(m_sampling_frequency * m_win_shift_ms / 1000);
}

void bob::ap::FrameExtractor::initWinSize()
{
  m_win_size = (size_t)pow(2.0,ceil(log((double)m_win_length)/log(2)));
  m_cache_frame_d.resize(m_win_size);
}

void bob::ap::FrameExtractor::initMaxRange()
{
  // update m_max_range, since m_sampling_frequency may have changed or set inside an Init()
  m_max_range = pow(2.0, m_sampling_frequency/1000)/2.0 - 0.5;
}

void bob::ap::FrameExtractor::extractNormalizeFrame(const blitz::Array<double,1>& input,
  const size_t i, blitz::Array<double,1>& frame_d) const
{
  // Set padded frame to zero
  frame_d = 0.;
  // Extract frame input vector
  blitz::Range rf(0,(int)m_win_length-1);
  blitz::Range ri(i*(int)m_win_shift,i*(int)m_win_shift+(int)m_win_length-1);
  frame_d(rf) = input(ri);

  if (m_normalize_mean) { // added by Pavel Korshunov
    // We normalize by subtracting mean value
    frame_d(rf) -= blitz::mean(frame_d);
  }
  else {
    //Otherwise, we normalize by dividing by maximum possible range, which is set in initWinLength()
    //This method of normalization is used in the following paper from Interspeech 2015:
    //"A Comparison of Features for Synthetic Speech Detection" by Md Sahidullah, Tomi Kinnunen, Cemal Hanilci
    if (m_max_range == 0)
      throw std::runtime_error("FrameExtractor: the maximum range in frame is 0. Please make sure you provide non-zero sampling frequency.");
    frame_d /= m_max_range;
  }
}


blitz::TinyVector<int,2>
bob::ap::FrameExtractor::getShape(const size_t input_size) const
{
  // Res will contain the number of frames x the dimension of the feature vector
  blitz::TinyVector<int,2> res;

  // 1. Number of frames
  res(0) = 1+((input_size-m_win_length)/m_win_shift);

  // 2. Dimension of the feature vector
  res(1) = m_win_length;

  return res;
}

blitz::TinyVector<int,2>
bob::ap::FrameExtractor::getShape(const blitz::Array<double,1>& input) const
{
  return getShape(input.extent(0));
}

