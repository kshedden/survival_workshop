# Survival analysis of NHANES III data
#
# Data sources:
# https://wwwn.cdc.gov/nchs/nhanes/nhanes3/datafiles.aspx
# https://www.cdc.gov/nchs/data-linkage/mortality-public.htm

import statsmodels.api as sm
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

pdf = PdfPages("nhanes.pdf")

dpath = "/nfs/kshedden/NHANES/"

# Read the survival data
fname = "NHANES_III_MORT_2011_PUBLIC.dat.gz"
colspecs = [(0, 5), (14, 15), (15, 16), (43, 46), (46, 49)]
names = ["seqn", "eligstat", "mortstat", "permth_int", "permth_exam"]
f = os.path.join(dpath, fname)
surv = pd.read_fwf(f, colspecs=colspecs, names=names, compression="gzip")

# Read the interview/examination data
fname = "adult.dat.gz"
colspecs = [(0, 5), (14, 15), (17, 19), (28, 31), (33, 34), (32, 33), (34, 35), (35, 41)]
names = ["seqn", "sex", "age", "county", "urbanrural", "state", "region", "poverty"]
f = os.path.join(dpath, fname)
df = pd.read_fwf(f, colspecs=colspecs, names=names, compression="gzip")
df = pd.merge(surv, df, left_on="seqn", right_on="seqn")

# Some data cleaning
df["poverty"] = df["poverty"].replace({888888: np.nan})
df["female"] = (df.sex == 2).astype(np.int)
df["rural"] = (df.urbanrural == 2).astype(np.int)
df["age_int"] = 12*df.age  # months
df["end"] = df.age_int + df.permth_int  # months

df = df.dropna()

# SurvfunRight can't handle 0 survival times
df = df.loc[df.end > df.age_int]


# The hazard function is the derivative of the cumlative hazard function.
def hazard(sf):
    tm = s.surv_times
    pr = s.surv_prob
    ii = (pr > 0)
    tm = tm[ii]
    pr = pr[ii]
    lpr = np.log(pr)
    return tm[0:-1], -np.diff(lpr) / np.diff(tm)

# Plot hazard functions for women and men
plt.clf()
plt.grid(True)
sex = {0: "Male", 1: "Female"}
for female in (0, 1):
    ii = df.female == female
    s = sm.SurvfuncRight(df.loc[ii, "end"], df.loc[ii, "mortstat"], entry=df.loc[ii, "age_int"])
    tm, hz = hazard(s)
    ha = sm.nonparametric.lowess(np.log(hz), tm/12)
    plt.plot(ha[:, 0], ha[:, 1], lw=3, label=sex[female])
ha, lb = plt.gca().get_legend_handles_labels()
leg = plt.figlegend(ha, lb, "upper center", ncol=2)
leg.draw_frame(False)
plt.xlabel("Age", size=15)
plt.ylabel("Log hazard", size=15)
plt.xlim(18, 90)
pdf.savefig()
pdf.close()

#
# Model 1
#

fml = "end ~ female + rural + C(region) + poverty"
model1 = sm.PHReg.from_formula(fml, status="mortstat", entry=df.age_int,
                               data=df)
result1 = model1.fit()

#
# Model 2
#

fml = "end ~ female + rural + C(region) + poverty"
model2 = sm.PHReg.from_formula(fml, status="mortstat", entry=df.age_int,
                               strata=df.state, data=df)
result2 = model2.fit()

#
# Model 3
#

fml = "end ~ female + rural + C(region) + poverty"
model3 = sm.PHReg.from_formula(fml, status="mortstat", entry=df.age_int,
                               strata=df.county, data=df)
result3 = model3.fit()
