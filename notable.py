# Survival analysis for the lifespans of notable people.  The data are
# available here:
#
# http://science.sciencemag.org/content/suppl/2014/07/30/345.6196.558.DC1

import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Run this once to get a simplified data set.
if False:
    df = pd.read_excel("SchichDataS1_FB.xlsx")
    df = df[["PrsLabel", "BYear", "DYear", "Gender"]]
    df.to_csv("schich.csv.gz", index=None, compression="gzip")

pdf = PdfPages("notable.pdf")

df = pd.read_csv("schich.csv.gz")
df = df.dropna()

df["lifespan"] = df.DYear - df.BYear

df = df.loc[df.Gender.isin(["Female", "Male"])]
df = df.loc[df.lifespan > 0]

#
# Survival function estimates for females and males
#

plt.clf()
ax = plt.axes([0.1, 0.1, 0.7, 0.8])
plt.grid(True)
plt.xlabel("Lifespan (years)", size=15)
plt.ylabel("Proportion", size=15)
for sex in "Female", "Male":
    ii = df.Gender == sex
    s = sm.SurvfuncRight(df.loc[ii, "lifespan"], np.ones(ii.sum()), title=sex)
    s.plot(ax=ax)

# Create a legend
ha, lb = ax.get_legend_handles_labels()
ha = [ha[0], ha[2]] # Optional, hide points from legend
lb = [lb[0], lb[2]]
leg = plt.figlegend(ha, lb, loc="center right")
leg.draw_frame(False)
pdf.savefig()

#
# Hazard estimates for females and males
#

# The hazard function is the derivative of the cumlative hazard function.
def hazard(sf):
    tm = s.surv_times
    pr = s.surv_prob
    ii = (pr > 0)
    tm = tm[ii]
    pr = pr[ii]
    lpr = np.log(pr)
    return tm[0:-1], -np.diff(lpr) / np.diff(tm)

plt.clf()
plt.grid(True)
for sex in "Female", "Male":
    ii = df.Gender == sex
    s = sm.SurvfuncRight(df.loc[ii, "lifespan"], np.ones(ii.sum()), title=sex)
    tm, hz = hazard(s)
    plt.plot(tm, np.log(hz), lw=3, label=sex)
ha, lb = plt.gca().get_legend_handles_labels()
leg = plt.figlegend(ha, lb, "upper center", ncol=2)
leg.draw_frame(False)
plt.xlabel("Age", size=15)
plt.ylabel("Log hazard", size=15)
plt.xlim(0, 90)
pdf.savefig()

#
# Survival function estimates by era
#

plt.clf()
ax = plt.axes([0.08, 0.1, 0.68, 0.8])
plt.grid(True)
plt.xlabel("Lifespan (years)", size=15)
plt.ylabel("Proportion", size=15)
for byear in np.arange(0, 2000, 500):
    ii = (df.BYear >= byear) & (df.BYear < byear + 500)
    s = sm.SurvfuncRight(df.loc[ii, "lifespan"], np.ones(ii.sum()),
                         title="%d-%d" % (byear, byear+500))
    s.plot(ax=ax)

# Create a legend
ha, lb = ax.get_legend_handles_labels()
ha = [ha[i] for i in range(0, len(ha), 2)] # Optional, hide points from legend
lb = [lb[i] for i in range(0, len(lb), 2)]
leg = plt.figlegend(ha, lb, loc="center right")
leg.draw_frame(False)
pdf.savefig()

#
# Hazard function estimates by era
#

plt.clf()
plt.grid(True)
for byear in np.arange(0, 2000, 500):
    ii = (df.BYear >= byear) & (df.BYear < byear + 500)
    s = sm.SurvfuncRight(df.loc[ii, "lifespan"], np.ones(ii.sum()))
    tm, hz = hazard(s)
    plt.plot(tm, np.log(hz), lw=3, label="%d-%d" % (byear, byear+500))
ha, lb = plt.gca().get_legend_handles_labels()
leg = plt.figlegend(ha, lb, "upper center", ncol=4)
leg.draw_frame(False)
plt.xlabel("Age", size=15)
plt.ylabel("Log hazard", size=15)
plt.xlim(0, 90)
pdf.savefig()

#
# Model 1: linear model for lifespan
#

# Question: when were the mean lifespans of women and men equal?
# In year y, women's mean lifespan was ... years different than
# mens' mean lifespan.
fml = "lifespan ~ BYear + Gender + BYear*Gender"
model1 = sm.OLS.from_formula(fml, data=df)
result1 = model1.fit()

for left in -1500, 500, 1500:
    plt.clf()
    plt.grid(True)
    plt.plot(df.BYear, df.lifespan, 'o', color='purple', alpha=0.5,
             rasterized=True)
    plt.xlabel("Year of birth", size=15)
    plt.ylabel("Lifespan", size=15)
    plt.xlim(left, 2000)
    pdf.savefig()

#
# Model 2: linear model for log lifespan
#

# In year y, women's mean lifepan was ... times greater/lesser than
# mens' mean lifespan.
fml = "I(np.log(lifespan)) ~ BYear + Gender + BYear*Gender"
model2 = sm.OLS.from_formula(fml, data=df)
result2 = model2.fit()


#
# Model 3: proportional hazards regression
#

fml = "lifespan ~ BYear + Gender + BYear*Gender"
model3 = sm.PHReg.from_formula(fml, data=df)
result3 = model3.fit()

# Plot the baseline hazard function
bhaz = result3.baseline_cumulative_hazard[0]
x = bhaz[0]
y = bhaz[1]
haz = np.diff(y, 1) / np.diff(x, 1)
plt.clf()
plt.grid(True)
plt.plot(x[0:-1], haz, lw=3)
plt.xlim(0, 90)
plt.ylim(0, 2)
plt.xlabel("Lifespan (years)", size=15)
plt.ylabel("Hazard", size=15)
pdf.savefig()

pdf.close()
