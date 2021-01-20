# Load the data in
data ./interval0pc.pi
# Ignore any "bad" channels
ignore bad
ignore **-0.3
ignore 10.0-**

# Using the model defined in the downloaded tarfile
@./models/interval0pc.xcm

# Apply the fit
fit

setplot rebin 100 50
iplot eeufspec

# The following have to run manually
hardcopy xrtSedPC.ps/ps
wd xrtSedPC.dat
