''' Present an interactive function explorer with slider widgets.
Scrub the sliders to change the properties of the ``sin`` curve, or
type into the title text box to update the title of the plot.
Use the ``bokeh serve`` command to run the example by executing:
    bokeh serve sliders.py
at your command prompt. Then navigate to the URL
    http://localhost:5006/sliders
in your browser.
'''
import numpy as np

from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Slider, TextInput, PreText, DataTable
from bokeh.plotting import figure, curdoc, output_file, save

from pytmatrix.tmatrix import Scatterer
from pytmatrix.psd import PSDIntegrator, GammaPSD
from pytmatrix import orientation, radar, tmatrix_aux, refractive

import pandas as pd
from scipy.special import gamma

# calculations
def drop_ar(D_eq):
    if D_eq < 0.7:
        return 1.0;
    elif D_eq < 1.5:
        return 1.173 - 0.5165*D_eq + 0.4698*D_eq**2 - 0.1317*D_eq**3 - \
            8.5e-3*D_eq**4
    else:
        return 1.065 - 6.25e-2*D_eq - 3.99e-3*D_eq**2 + 7.66e-4*D_eq**3 - \
            4.095e-5*D_eq**4 

def get_scattering_props(Dm,logNw,mu):
    scatterer = Scatterer(wavelength=tmatrix_aux.wl_S, m=refractive.m_w_10C[tmatrix_aux.wl_S])
    scatterer.psd_integrator = PSDIntegrator()
    scatterer.psd_integrator.axis_ratio_func = lambda D: 1.0/drop_ar(D)
    scatterer.psd_integrator.D_max = 10.0
    scatterer.psd_integrator.geometries = (tmatrix_aux.geom_horiz_back, tmatrix_aux.geom_horiz_forw)
    scatterer.or_pdf = orientation.gaussian_pdf(20.0)
    scatterer.orient = orientation.orient_averaged_fixed
    scatterer.psd_integrator.init_scatter_table(scatterer)

    D0 = (3.67 + mu)/(4 + mu) * Dm
    Nw = 10.**logNw   
    scatterer.psd = GammaPSD(D0=D0, Nw=Nw, mu=mu)

    Zh = 10*np.log10(radar.refl(scatterer))
    Zv = 10*np.log10(radar.refl(scatterer, False))
    Zdr = radar.Zdr(scatterer)
    Ldr = 10*np.log10(radar.ldr(scatterer))
    rho_hv = radar.rho_hv(scatterer)
    scatterer.set_geometry(tmatrix_aux.geom_horiz_forw)
    Kdp = radar.Kdp(scatterer)
    Ah = radar.Ai(scatterer)



    f_u = (6/(4**4))*((4 + mu)**(mu + 4))/(gamma(mu+4))

    dsd_df = pd.DataFrame()
    dsd_df['D'] = np.arange(0.1,10,0.1)
    dsd_df['ND'] = Nw * f_u * ((dsd_df['D'] / Dm) ** mu) * np.exp(-1.*(4 + mu)*(dsd_df['D'] / Dm))
    NT = Nw * f_u * gamma(mu + 1) * Dm / ((4 + mu)**(mu + 1))
    LWC = (1. * np.pi * Nw * Dm**4) / (4.**4 * 1000.)

    integ_df = pd.DataFrame({'NT':[NT], 
                            'LWC':[LWC], 
                            'Zh':[Zh], 
                            'Zv':[Zv], 
                            'Zdr': [Zdr], 
                            'Ldr': [Ldr], 
                            'rho_hv':[rho_hv], 
                            'Kdp':[Kdp], 
                            'Ah':[Ah]})

    return dsd_df, integ_df

def calc_dsd(Dm,logNw,mu):
    f_u = (6/4**4)*((4 + mu)**(mu + 4))/(gamma(mu+4))
    dsd_df = pd.DataFrame()
    dsd_df['D'] = np.arange(0.1,10,0.1)
    dsd_df['ND'] = 10**logNw * f_u * ((dsd_df['D'] / Dm) ** mu) * np.exp(-1.*(4 + mu)*(dsd_df['D'] / Dm))
 
    return dsd_df

def update_data(attrname, old, new):

    # Get the current slider values
    Dm_new = Dm_slider.value
    logNw_new = logNw_slider.value
    mu_new = mu_slider.value

    # Generate the new curve
    dsd_df_new, integ_df_new = get_scattering_props(Dm_new, logNw_new, mu_new)
#    x_new = dsd_df_new['D']
#    y_new = dsd_df_new['ND']
#    print(y_new)
    update_stats(integ_df_new)

    source.data = dsd_df_new
    source2.data = integ_df_new

def update_stats(integ_df):
    stats.text = str(integ_df.T)


# Set up data
Dm = 2.0
logNw = 3.0
mu = 0
dsd_df, integ_df = get_scattering_props(Dm, logNw, mu)
#x = dsd_df['D']
#y = dsd_df['ND']

source2 = ColumnDataSource(data=integ_df)

stats = PreText(text='', width=500)
#stats = DataTable(source=source)

source = ColumnDataSource(data=dsd_df)

update_stats(integ_df)

# Set up plot
plot = figure(height=400, width=400, title="Particle size distribution",
              tools="crosshair,pan,reset,save,wheel_zoom",y_axis_type="log",
              x_range=[0, 10], y_range=[1,1000000])

plot.line('D', 'ND', source=source, line_width=3, line_alpha=0.6)


# Set up widgets
text = TextInput(title="title", value='my PSD')
Dm_slider = Slider(title="Dm", value=3.0, start=0, end=6.0, step=0.5)
logNw_slider = Slider(title="log_10 Nw", value=3, start=0.5, end=6, step=0.5)
mu_slider = Slider(title="mu", value=0.0, start=-3.5, end=20., step=0.5)


for w in [Dm_slider, logNw_slider, mu_slider]:
    w.on_change('value_throttled', update_data)


# Set up layouts and add to document
inputs = column(text, Dm_slider, logNw_slider, mu_slider)

curdoc.add_root(row(inputs, plot, stats, width=800))
curdoc.title = "My PSD Explorer"

curdoc.theme = Theme(json=yaml.load("""
    attrs:
        Figure:
            background_fill_color: "#DDDDDD"
            outline_line_color: white
            toolbar_location: above
            height: 500
            width: 800
        Grid:
            grid_line_dash: [6, 4]
            grid_line_color: white
""", Loader=yaml.FullLoader))

