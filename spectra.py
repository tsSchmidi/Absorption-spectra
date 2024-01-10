import pandas as pd
import pylab
import scipy.optimize as opt
from scipy.interpolate import interp1d

################################## input #####################################

data = pd.read_excel(r"C:\Users\s1451723\Desktop\Stuff\Spectra\Raw data\test3.xlsx")
figure_rows = 1 #Output figure dimenions
figure_columns = 1
noempty = True #Don't plot empty wells
title = ""
ylab = "delta normalised OD"
plotsample = False #List of sample and replicate numbers, False if want all
plotreplicate = False
plotmean = True
plotcontrol = True
commonscale = False
mode = "line norm"
control = "EK01"
blank = "Water"
outliers = []
wl_min = 400
wl_max = 660
wl_norm = 450

################################ Functions ####################################

def func(x, a, b):
    return(a*x + b)

############################### Prepare data ##################################

#Trim data
data = data.drop([wl for wl in list(data)[1:] if wl < wl_min or wl > wl_max], axis=1)
samples = list(dict.fromkeys([data.values[i][0] for i in range(len(data))]))
samples = [sample for sample in samples if sample != blank]
if plotsample:
    samples = [samples[i-1] for i in plotsample]
if blank:
    blank_curve = data.groupby([list(data)[0]]).get_group(blank).mean()
data = data[data[list(data)[0]].isin(samples+[control])]

for sample in outliers:
    sample_mean = data.groupby([list(data)[0]]).get_group(sample).mean()
    d = []
    replicates = data[data[list(data)[0]] == sample]

    for replicate in [replicates[list(replicates)[1:]].iloc[i] for i in range(len(replicates))]:
        d += (sample_mean - replicate).abs().max(),
    data = data.drop(index = replicates.index[[i == max(d) for i in d].index(True)])

#Split sample names from data
data1 = data[list(data)[0]]
data2 = data[list(data)[1:]].copy()

#Subtract blank
if blank:
    data2 = data2 - blank_curve

#######################
if mode == "line norm":
    
    fitx = list(data2)[:3]+list(data2)[-3:]

    for i in range(len(data2)):
        fity = list(data2.iloc[i])[:3]+list(data2.iloc[i])[-3:]
        a,b = list(opt.curve_fit(func,fitx,fity)[0])
        data2.iloc[i] += [-func(x,a,b) for x in list(data2)]
        #data2.iloc[i] += -min(data2.iloc[i])
        
######

if mode == "OD norm":
    #Divide everything with OD
    data2 = data2.div(data2[wl_norm],axis=0)

######

if mode == "interpolation":
    # requires multiple concentrations of control
    x = list(data2)
    data_control = data[data[list(data)[0]]==control][list(data)[1:]]
    interps = [interp1d(
        data_control[wl_norm],
        data_control[wl],
        fill_value="extrapolate")
               for wl in x]
    interps = [interp(data2[wl_norm]) for interp in interps]
    interps = pd.DataFrame(interps).transpose()
    interps.columns = data2.columns
    interps.index = data2.index
    data2 = data2.sub(interps)
#######################
    
#Rejoin data
data = pd.concat([data1,data2],axis=1)

#######################
if mode != "raw":
    #Construct a control curve based on WT (need sample names for this)
    control_curve = data.groupby([list(data)[0]]).get_group(control).mean()

    #Split again
    data1 = data[list(data)[0]]
    data2 = data[list(data)[1:]]

    #Subtract control from everything
    data2 = data2 - control_curve

    #Rejoin again
    data = pd.concat([data1,data2],axis=1)
#######################

#Get min and max values, can be used to set ylim
mini = min([min(i) for i in data2.values])
maxi = max([max(i) for i in data2.values])

ranges = data.groupby([list(data)[0]]).max().max(axis=1).subtract(
    data.groupby([list(data)[0]]).min().min(axis=1))
if ranges[control] != ranges.min():
    minrange = ranges[control]
else:
    minrange = ranges.drop(control).min()

################################ Plot data ####################################

#x: wavelengths, z: unique sample names
x = list(data)[1:]

#Get x label from excel, y label manually
xlab = list(data)[0]

#Empty canvas in inches
fig = pylab.figure(figsize = (min(12, 6*figure_columns/figure_rows),min(6,12*figure_rows/figure_columns)))

#Shared X and Y axis labels
fig.text(0.5, 0.02, xlab, ha='center')
fig.text(0.01, 0.5, ylab, va='center', rotation='vertical')
fig.text(0.5, 0.95, title, ha='center', fontsize=12)

#Count plots for position
n = 1

for sample in samples:
    #Dictates the destination of the next plot
    fig.add_subplot(figure_rows, figure_columns, n)
    n += 1

    #Get and loop through set of rows that have a given sample name
    replicates = data[data[list(data)[0]] == sample]
    ind = list(range(len(replicates)))
    if plotreplicate:
        ind = [ind[i-1] for i in plotreplicate]

    if plotmean:
        y = replicates.mean().values
        #bar = [replicates.mean().subtract(replicates.min(numeric_only=True)).values, replicates.max(numeric_only=True).subtract(replicates.mean()).values]
        #pylab.errorbar(x,y,bar)
        pylab.plot(x,y)

    for replicate in [replicates[i:i+1] for i in ind]:
        #y: delta normalised OD compared to control
        y = replicate.drop(list(data)[0], axis=1).values[0]
        pylab.scatter(x,y,s=3)
        #End of replicate

    if mode != "raw" and plotcontrol:
        #Straight line y=0
        pylab.plot(x,[0 for i in x], label="Control")
        pylab.legend(loc="upper right")

    #Subplot title and limits
    pylab.title(sample)
    if commonscale:
        pylab.ylim(mini-(maxi-mini)/40,maxi+(maxi-mini)/40)
    elif ranges[sample] < minrange:
        pad = minrange - ranges[sample]
        ymin, ymax = pylab.ylim()
        pylab.ylim(ymin - pad/2, ymax + pad/2)
    #End of sample

#Prevent overlapping
pylab.tight_layout(pad=2)
pylab.show()
