/**
 * @author Elie Khoury <Elie.Khoury@idiap.ch>
 * @author Laurent El Shafey <Laurent.El-Shafey@idiap.ch>
 * @date Wed Jan 11:10:20 2013 +0200
 *
 * @brief Implement a rectangular window frame extractor
 *
 * Copyright (C) Research Institute, Martigny, Switzerland
 */

#ifndef BOB_AP_FRAME_EXTRACTOR_H
#define BOB_AP_FRAME_EXTRACTOR_H

#include <blitz/array.h>

namespace bob {
/**
 * \ingroup libap_api
 * @{
 *
 */
namespace ap {

/**
 * @brief This class is a base class for classes that perform audio processing
 * on a frame basis.
 */
class FrameExtractor
{
  public:
    /**
     * @brief Constructor. Initializes working arrays
     */
    FrameExtractor(const double sampling_frequency,
      const double win_length_ms=20., const double win_shift_ms=10.,
      const bool normalize_mean=true);

    /**
     * @brief Copy Constructor
     */
    FrameExtractor(const FrameExtractor& other);

    /**
     * @brief Assignment operator
     */
    FrameExtractor& operator=(const FrameExtractor& other);

    /**
     * @brief Equal to
     */
    bool operator==(const FrameExtractor& other) const;

    /**
     * @brief Not equal to
     */
    bool operator!=(const FrameExtractor& other) const;

    /**
     * @brief Destructor
     */
    virtual ~FrameExtractor();

    /**
     * @brief Gets the output shape for a given input/input length
     */
    virtual blitz::TinyVector<int,2> getShape(const size_t input_length) const;
    virtual blitz::TinyVector<int,2> getShape(const blitz::Array<double,1>& input) const;

    /**
     * @brief Returns the sampling frequency/frequency rate
     */
    double getSamplingFrequency() const
    { return m_sampling_frequency; }
    /**
     * @brief Returns the window length in miliseconds
     */
    double getWinLengthMs() const
    { return m_win_length_ms; }
    /**
     * @brief Returns the window length in number of samples
     */
    size_t getWinLength() const
    { return m_win_length; }
    /**
     * @brief Returns the window shift in miliseconds
     */
    double getWinShiftMs() const
    { return m_win_shift_ms; }
    /**
     * @brief Returns the window shift in number of samples
     */
    size_t getWinShift() const
    { return m_win_shift; }
    /**
     * @brief Tells whether frame should be normalized by subtracting mean (True) or dividing by max_range (False)
     */
    bool getNormalizeMean() const
    { return m_normalize_mean; }

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
     * @brief Sets whether frame should be normalized by subtracting mean (True) or dividing by max_range (False)
     */
    virtual void setNormalizeMean(const double normalize_mean)
    { m_normalize_mean = normalize_mean; }

  protected:
    /**
     * @brief Extracts the frame of the given index
     * @warning No check is performed
     */
    virtual void extractNormalizeFrame(const blitz::Array<double,1>& input,
      const size_t i, blitz::Array<double,1>& frame) const;
    virtual void initWinSize();
    virtual void initWinLength();
    virtual void initWinShift();
    virtual void initMaxRange();

    double m_sampling_frequency; ///< The sampling frequency
    double m_win_length_ms; ///< The window length in miliseconds
    size_t m_win_length;
    double m_win_shift_ms;
    size_t m_win_shift;
    size_t m_win_size;
    double m_max_range; //half of the maximum possible dynamic range of the original signal (for 16 bits, it is 32768)
    bool m_normalize_mean; //normalize the frame by subtracting its mean

    mutable blitz::Array<double,1> m_cache_frame_d;
};

}
}

#endif /* BOB_AP_FRAME_EXTRACTOR_H */
