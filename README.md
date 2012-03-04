Guitar Tuner
============

About
-----

This is my simple project to learn about Fourier Transforms. The Fourier
Transform converts a function from the time domain to the frequency
domain. This makes it easy to determine the dominant frequency in a
signal. A tuner only has to compare this frequency against known
frequencies for musical notes.

Building
--------

Requires:

   * alsa
   * fftw3
   * sdl (optional, used for displaying the frequency histogram)

Eventually, I plan to write my own implementation of a Fast Fourier
Transform algorithm once I understand it. Currently I have the basic
implementation of the Fourier Transform implemented, but it is generally
too slow to use as it is O(n^2).


