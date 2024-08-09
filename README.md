# cdrgam
An `R` suite of utility functions to support continuous-time
deconvolutional regression (CDR) of nonlinear impulse response
functions (IRFs), based on generalized additive models (GAMs) as
implemented by Simon Woods' `mgcv` package. No modifications are made
to the model itself, since `mgcv`'s native "linear functional terms"
feature can be exploited to ensure that a GAM smooth implements a
continuous-time and potentially nonlinear convolution kernel. In other
words, a CDR model can be instantiated as a GAM model using
functionality already available through `mgcv`. However, in practice,
implementation is non-trivial: in order for a GAM to be a valid CDR
model, the data must be organized in a particular way and the model
must have a particular structure.  Moreover, `mgcv`'s standard
visualization tools are ill-suited to CDR applications. This package
facilitates the use of `mgcv` for CDR modeling by providing helper
functions to organize data, define models, visualize estimates, and
perform statistical model comparison. These functions reduce the
likelihood of user error and simplify the process of spinning up new
CDR-GAM analyses from scratch.
