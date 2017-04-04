Survival Analysis in Python Statsmodels
=======================================

Source code for statsmodels survival methods:

https://github.com/statsmodels/statsmodels/tree/master/statsmodels/duration

Survival analysis is used to analyze data in which the primary
variable of interest is a time duration.  For example, the time
duration could be: a person's lifespan, the time that a person
survives after being diagnosed with a serious disease, the time from
the diagnosis of a disease until the person recovers, or the duration
of time that a piece of machinery remains in good working order.

Duration data can be used to answer questions analogous to what we
might consider for other types of data: what is the mean duration of a
population?; what is the 75th percentile of the durations in a
population?; when comparing two populations, which one has the shorter
expected or median duration?; can we identify the unique associations
between many independent variables and a duration outcome?

Here are some key concepts in duration analysis:

* __Time origin__: The durations of interest correspond to time
  intervals that begin and end when something specific happens.  It is
  usually important to be very explicit about what defines the time
  origin (time zero) from which the duration is calculated.  For
  example, when looking at the survival of people with a disease, the
  time origin could be the date of diagnosis.  When looking at human
  life spans ("all cause mortality") it might make more sense to
  define the time origin to be the date of birth.

* __Event__: This term refers to whatever happens that concludes the
  time interval of interest.  It may be death, or some other type of
  "failure", or it may be something more favorable, like recovery or
  cure.  Most survival analysis is based on the idea that every
  subject will eventually experience the event.

* __Survival time distribution__: This is a marginal distribution
  defining the proportion of the population that has experienced the
  event on or before time T.  Usually it is expressed as the
  complementary "survival function" (e.g. the proportion of people who
  have not yet died).

* __Censoring__: Censoring occurs when we do not observe when a
  subject experiences the event of interest, but we do have some
  partial information about that time.  The most common form of
  censoring is _right censoring_, in which we observe a time T such
  that we know the event did not occur prior to time T.

* __Risk set__: This is the set of units (e.g. people) in a sample at
  a given time who may possibly experience the event at that time.  It
  is usually the set of people who have not already experienced the
  event or been censored (but may be only a subset of these people).

* __Hazard__: This is the probability of experiencing the event in the
  next time unit, given that it has not already occurred (technically,
  this is the discrete time definition, the continuous time definition
  involves rates but follows the same logic).


Survival function estimation
----------------------------

If there is no censoring, it is extremely easy to estimate the
survival function.  It is simply the complement of the empirical
cumulative distribution function of the data.

If there is censoring, the standard method for estimating the survival
function is the product-limit estimator or [Kaplan Meier
estimator](https://en.wikipedia.org/wiki/Kaplan%E2%80%93Meier_estimator).

The idea behind the Kaplan-Meier estimate is not difficult.  Group the
data by the distinct times t(1) < t(2) < ..., and let R(t) denote the
risk set size at time t, and let d(t) indicate the number of events at
time t (if there are no ties, d(t) will be equal to either 0 or 1).
The probability of the event occuring at time t is estimated to be
d(t)/R(t).  The probability of the event not occurring at time t is
therefore estimated to be 1 - d(t)/R(t).  The probability of making it
to time t wihout experiencing the event is therefore estimated to be

(1 - d(1)/R(1)) * (1 - d(2)/R(2)) * ... * (1 - d(t)/R(t)).


Regression analysis for survival/duration data
----------------------------------------------

Regression analysis is used to understand how multiple factors of
interest are related to an "outcome" variable of interest.  If this
outcome variable is a duration, we are doing "survival regression".

By direct analogy with linear regression, we might seek to model the
expected survival time as a function of covariates.  If there is no
censoring, we could, for example, use least square regression to
relate the survival time T, or a transformation of it (e.g. log(T)) to
a linear function of the covariates.  While this is sometimes done, it
is more common to approach regression for duration data by modeling
the hazard rather than modeling the duration itself.

The hazard is the conditional probability of experiencing the event of
interest at time T, given that it has not yet occurred.  For example,
in a medical study this may be the probability of a subject dying at
time T given that the subject was still alive just before time T (in
continuous time we would substitute "rate" for "probability" but will
ignore this distinction here).

In survival regression, we view the hazard as a function that is
determined by the covariates.  For example, the hazard may be
determined by age and gender.  A very popular form of hazard
regression models the conditional hazards as a collection of parallel
functions, specifically

h(t, x) = b(t) * exp(c0 + c1*x1 + ... + cp*xp),

where b(t) is the "baseline hazard function", the scalars c0, ..., cp
are unknown regression coefficients, and the xj are the observed
covariates for one subject.  This model can also be written in log
form

log h(t, x) = log b(t) + c'x

where c=(c0, ..., cp), x=(1, x1, ..., xp).  Thus, the log hazard is
modeled as a time varying intercept plus a linear predictor that is
not time varying (there are generalizations in which the linear
predictor is also time varying).

This regression model is called "proportional hazards regression" or
the "Cox model".  A key feature of this model is that it is possible
to estimate the coefficients c using a partial likelihood that does
not involve the baseline hazard function.  This makes the procedure
"semi-parametric".

The key point to remember about interpreting this model is that a
coefficient cj, for a covariate, say age, has the property that the
hazard (e.g. of dying) changes multiplicatively by a factor of exp(cj)
for each unit increase in the covariate's value.
