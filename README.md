Survival Analysis in Python Statsmodels
=======================================

Source code for statsmodels survival methods:

https://github.com/statsmodels/statsmodels/tree/master/statsmodels/duration

Introduction
------------

[Survival analysis](https://en.wikipedia.org/wiki/Survival_analysis)
is used to analyze data in which the primary variable of interest is a
time duration.  For example, the time duration could be a person's
lifespan, the time that a person survives after being diagnosed with a
serious disease, the time from the diagnosis of a disease until the
person recovers, or the duration of time that a piece of machinery
remains in good working order.

Duration data can be used to answer questions analogous to what we
might consider for other types of data, such as: what is the mean
duration of a population?; what is the 75th percentile of the
durations in a population?; when comparing two populations, which one
has the shorter expected or median duration?; can we identify the
unique associations between many independent variables and a duration
outcome?

Here are some key concepts in duration analysis:

* __Time origin__: The durations of interest correspond to time
  intervals that begin and end when something specific happens.  It is
  important to be very explicit about what defines the time origin
  (time zero) from which the duration is calculated.  For example,
  when looking at the survival of people with a disease, the time
  origin could be the date of diagnosis.  When looking at human
  lifespans ("all cause mortality") it might make more sense to define
  the time origin to be the date of birth.

* __Event__: This term refers to whatever happens that concludes the
  time interval of interest.  It may be death, or some other type of
  "failure", or it may be something more favorable, like recovery or
  cure from a disease.  Most survival analysis is based on the idea
  that every subject will eventually experience the event.

* __Survival time distribution__: This is a marginal distribution
  defining the proportion of the population that has experienced the
  event on or before time T.  Usually it is expressed as the
  complementary "survival function" (e.g. the proportion of people who
  have not yet died as of time T).

* __Censoring__: Censoring occurs when we do not observe when a
  subject experiences the event of interest, but we do have some
  partial information about that time.  The most common form of
  censoring is _right censoring_, in which we observe a time T such
  that we know the event did not occur prior to time T.  Other forms
  of censoring are _interval censoring_ and _left censoring_.

* __Risk set__: This is the set of units (e.g. people) in a sample at
  a given time who may possibly experience the event at that time.  It
  is usually the set of people who have not already experienced the
  event and have not been censored (but the risk set may be only a
  subset of these people).

* __Hazard__: This is the probability of experiencing the event in the
  next time unit, given that it has not already occurred (technically,
  this is the discrete time definition, the continuous time definition
  of the hazard involves rates but follows the same logic).


Marginal survival function and hazard estimation
------------------------------------------------

If there is no censoring, it is extremely easy to estimate the
marginal survival function.  It is simply the complement of the
empirical [cumulative distribution
function](https://en.wikipedia.org/wiki/Cumulative_distribution_function)
of the data.

If there is censoring, the standard method for estimating the survival
function is the product-limit estimator or [Kaplan Meier
estimator](https://en.wikipedia.org/wiki/Kaplan%E2%80%93Meier_estimator).

The idea behind the Kaplan-Meier estimate is not difficult.  Group the
data by the distinct times t(1) < t(2) < ..., and let R(t) denote the
risk set size at time t, and let d(t) indicate the number of events at
time t (if there are no ties, d(t) will always be equal to either 0 or
1).  The probability of the event occurring at time t (given that it
has not occurred already) is estimated to be d(t)/R(t).  The
probability of the event not occurring at time t is therefore
estimated to be 1 - d(t)/R(t).  The probability of making it to time t
without experiencing the event is therefore estimated to be

(1 - d(1)/R(1)) ⋅ (1 - d(2)/R(2)) ⋅ ... ⋅ (1 - d(t)/R(t)).

A consequence of this definition is that the estimated survival
function obtained using the product limit method is a step function,
with steps at the event times.

A closely related
[calculation](https://en.wikipedia.org/wiki/Nelson%E2%80%93Aalen_estimator)
estimates the marginal hazard function.

Survival and hazard functions are usually presented as plots

Regression analysis for survival/duration data
----------------------------------------------

Regression analysis is used to understand how multiple factors of
interest are related to an "outcome" variable of interest.  If this
outcome variable is a duration, we are doing "survival regression".

By direct analogy with linear regression, we might seek to model the
expected survival time as a function of covariates.  If there is no
censoring, we could, for example, use least squares regression to
relate the survival time T, or a transformation of it (e.g. log(T)) to
a linear function of the covariates.  While this is sometimes done, it
is more common to approach regression for duration data by modeling
the hazard rather than modeling the duration itself.

The hazard is the conditional probability of experiencing the event of
interest at time T, given that it has not yet occurred.  For example,
in a medical study this may be the probability of a subject dying at
time T given that the subject was still alive at time T (in continuous
time we would substitute "rate" for "probability" but we ignore this
distinction here).

In survival regression, we view the hazard as a function that is
determined by the covariates.  For example, the hazard may be
determined by age and gender.  A very popular form of hazard
regression models the conditional hazards as a collection of parallel
functions, specifically

h(t, x) = b(t) ⋅ exp(c0 + c1⋅x1 + ... + cp⋅xp),

where b(t) is the "baseline hazard function", the scalars c0, ..., cp
are unknown regression coefficients, and the xj are the observed
covariates for one subject.  This model can also be written in log
form

log h(t, x) = log b(t) + c'x

where c = [c0, ..., cp], x = [1, x1, ..., xp].  Thus, the log hazard
is modeled as a time varying intercept plus a linear predictor that is
not time varying (there are generalizations in which the linear
predictor is also time varying).

This regression model is called [proportional hazards
regression](https://en.wikipedia.org/wiki/Proportional_hazards_model)
or the "Cox model".  A key feature of this model is that it is
possible to estimate the coefficients c using a partial likelihood
that does not involve the baseline hazard function.  This makes the
procedure "semi-parametric".

The key point to remember about interpreting this model is that a
coefficient cj, for a covariate, say age, has the property that the
hazard (e.g. of dying) changes multiplicatively by a factor of exp(cj)
for each unit increase in the covariate's value.

More on censoring
-----------------

A large part of survival analysis is concerned with appropriately
handling censoring (if there is no censoring, it is generally possible
to analyze the log durations using standard, non-survival methods).
Censoring can be a subtle topic.  All survival methods have
limitations on the type of censoring they can handle, and it is not
always easy or even possible to determine in a given setting whether a
survival method can accommodate the type of censoring that is present.

To make things more concrete, we usually imagine that every subject
has both an event time T and a censoring time C.  That is, every
subject would eventually experience the event (if there were no
censoring), and would eventually be censored (if the event did not
happen).  We observe min(T, C), and an indicator of whether death
occured (i.e. that min(T, C) is equal to T).

The key requirement for most surival methods is that we have
"independent censoring", meaning that T and C are statistically
independent quantities.  Since we never observe both T and C for the
same person, it is usually not possible to directly assess whether T
and C are dependent.  However knowledge about the data collection
process is often used to assess whether independent censoring is
plausible to hold.

For example, one type of censoring that is quite common is
"administrative censoring".  This occurs when a study has a fixed data
collection window, say a thee year interval from January 1, 2012 to
January 1, 2015.  Suppose that people are randomly recruited into the
study, and if the event has not occured by January 1, 2015, the
subject is censored.  Thus, subjects who are recruited into the study
later are more likely to be censored.  But as long as the reruitment
date is not dependent with the true survival time, then T and C are
independent.  We can imagine a setting when administrative censoring
may induce dependence, e.g. if the subjects recruited later in the
study were healthier than those recruited earlier.  But in many
situations, this can be excluded as a likely circumstance based on
knowledge of how the study was conducted.

On the other hand, in some cases there is strong reason to believe
that subjects are more likely to be censored as they grow sicker,
which may mean that T and C are positively dependent.  For example, if
we have medical study in which the data come from insurance records,
as people get sicker they are more likely to become unable to work,
and may have to quit their job (leading to them being censored).
Similarly, as people age they may retire or become eligible for
Medicare, leading to age-dependent censoring.  Since age is likely
correlated with survival time, this could also induce dependent
censoring.

One remedy for dependent censoring is to identify covariates such that
T and C become independent after conditioning on the covariates.  For
example, age or a measure of overall health may be sufficient to
substantially reduce the dependence between T and C.