# Load the data in
data ./interval0wt.pi
# Ignore any "bad" channels
ignore bad
ignore **-0.3
ignore 10.0-**

# Using the model defined in the downloaded tarfile
@./models/interval0wt.xcm

# Apply the fit
fit

setplot rebin 100 50
iplot eeufspec

# The following have to run manually
hardcopy xrtSedWT.ps/ps
wd xrtSedWT.dat
