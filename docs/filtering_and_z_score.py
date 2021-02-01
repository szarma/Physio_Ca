#' % Event detection and processing
#' % Srdjan Sarikas
#'
#+echo=False
import numpy as np
import matplotlib.pyplot as plt
#' # Noise in microscopy
#'
#'  ## Ideal detector
#' Light comes packed in photons. In an ideal detector, the number of photons detected is only a linear function of the dwell time.
#' However, measuring exactly same light source twice, for the exact same time, does not need to yield the same number of photons detected, even for an ideal idetector.
#' A distribution of photon numbers obtained from many such measurements is called a [Poisson distribution](https://en.wikipedia.org/wiki/Poisson_distribution).
#' One of the properties of the Poisson distribution is that its standard deviation is the square root of its mean.
#' Suppose the source is shooting 0.01 photon per second on average
#'
#'
#'
#'
#'
#'
#'