#pragma once

#include <chrono>
#include <iostream>
#include <cmath>

namespace cbica
{
  /**
  \brief A class to handle progress bar in CLI applications

  Usage example:
  \verbatim
  float limit = 10000;
  cbica::ProgressBar progressBar(limit); // width of bar defaults to 100
  for(float i = 0; i < limit; i++)
  {
    ++progressBar; // increase counter
    progressBar.display(); // show the bar
    // do not have any other message coming here, otherwise the progress bar will get wiped
  }
  progressBar.done(); // wrap-up
  \endverbatim
  */
  class ProgressBar
  {
  private:
    float ticks = 0;

    const float total_ticks;
    const float bar_width;
    const char complete_char = '=';
    const char incomplete_char = ' ';
    const std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  public:
    //! Full constructor
    ProgressBar(float total, float width, char complete, char incomplete) :
      total_ticks{ total }, bar_width{ width }, complete_char{ complete }, incomplete_char{ incomplete } {}

    //! Partial constructor
    ProgressBar(float total, float width) : total_ticks{ total }, bar_width{ width } {}

    //! Partial constructor
    ProgressBar(float total) : total_ticks{ total }, bar_width{ 100 } {}

    unsigned int operator++() { return ++ticks; }

    void display() const
    {
      auto progress = ticks / total_ticks;
      int pos = std::round(bar_width * progress);

      std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
      auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time).count();

      std::cout << "[";

      for (int i = 0; i < bar_width; ++i) 
      {
        if (i < pos) std::cout << complete_char;
        else if (i == pos) std::cout << ">";
        else std::cout << incomplete_char;
      }
      std::cout << "] " << int(progress * 100.0) << "% "
        << float(time_elapsed) / 1000.0 << "s\r";
      std::cout.flush();
    }

    void done() const
    {
      display();
      std::cout << std::endl;
    }
  };
}
