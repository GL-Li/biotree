#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 22:28:40 2018

@author: gl
"""

def fit_backforce(df=force, sample="BT1567_SF"):
    # match to BT1567 mannually
#    plt.plot(brain.index, brain["BT1567_SF"], color = "orange", lw=0.5)
#    plt.ylim(0, 20)
#    plt.xlim(0, 5100)

    sect_0 = df.loc[1000:1430, [sample]]
    sect_1 = df.loc[1460:1630, [sample]]
    sect_2 = df.loc[1700:4500, [sample]]

    slope_0 = get_slope(sect_0)
    slope_1 = get_slope(sect_1)
    slope_2 = get_slope(sect_2)

    # section 0 fit
    x0 = np.arange(500, 1700)
    y0 = slope_0.iloc[0]["intercept"] + slope_0.iloc[0]["slope"] * x0
    #plt.plot(x0, y0, color = "red", ls="--")
    
    # section 1 fit
    x1 = np.arange(1000, 1750)
    y1 = slope_1.iloc[0]["intercept"] + slope_1.iloc[0]["slope"] * x1
    #plt.plot(x1, y1, color = "red", ls="--")
    
    # section 2 fit
    x2 = np.arange(1400, 5000)
    y2 = slope_2.iloc[0]["intercept"] + slope_2.iloc[0]["slope"] * x2
    #plt.plot(x2, y2, color = "red", ls="--")
           
    x = np.concatenate((x0, x1, x2))
    y = np.concatenate((y0, y1, y2))
#    plt.scatter(x, y, s=0.5, c="red")
    return(pd.DataFrame({"x" : x, "y" : y}))

test = fit_backforce()

# average density of each sample
def get_density():
    brain_vaso = convert_vaso_to_df("../brain_perfusion/vasometrics_data/all_fast_perfusion.csv", skip_rows=3)
    brain_vaso.columns = [w[0:6] for w in brain_vaso.columns]
    return(get_avg_density(brain_vaso))
 

def compare_to_BT1567(df, sample, x_max=6000, y_max=20, title="CD1 mouse", flow_rate=100, friction=0, harden_time=4750,
                      xshift=0, yshift=0, show_density=False):
    """
    plot the fitted lines of BT1567 to sample in df
    
    Parameters
    ----------
    df : dataframe containing the data
    sample: string
        name of the sample, for example "BT1560_SF"
    flow_rate: number
        flow rate in perfusion, in ml/hr
    friction: number
        friction of syringe, in N
    harden_time: number
        time data point when back force reaches 20 N
    xshif, yshit : number
        move the fitted lines in x and y direction
        
    Returns
    -------
    a figure of the sample's force curve and fitted lines of BT1567_SF
    """
    # adjust time based on harden time
    k_x = harden_time / 4750   # 4750 is for BT1567, or for 500 mg catalyst at 24 degree
    
    # flow rate compared to green syringe used for BT1567, used to adjust force
    k_y = flow_rate / 100
    
    # plot the adjusted experimental data according to parameters of BT1567
    # syringe friction is removed so the force is pure hydrodynamic resistence
    #plt.figure(figsize=(6, 4))
    plt.plot(df.index, (df[sample] - friction), color = "blue")
    
    # get fitting data to BT1567
    xy = fit_backforce()
    x = xy.x
    y = xy.y - 0.177  # subtract friction to get real hydro resistance
    plt.scatter(x * k_x + xshift, y * k_y + yshift, c="red", s=0.5)
    plt.ylim(-0.3, y_max)
    plt.xlim(0, x_max)
    
    plt.annotate(sample, xy=(100, y_max - 1.5), color="blue")
    # display density if available
    if show_density == True:
            # density of the sample
            density = round(densities[sample[0:6]])
            plt.annotate("density: " + str(density), xy=(100, 13), color="red")
    plt.annotate("x_shift: " + str(xshift), xy=(x_max-2000, 2))
    plt.annotate("y_shift: " + str(yshift), xy=(x_max-2000, 0.5))
    plt.title(title)