import astroquery
import numpy as np
import re
import os
from astroquery.alma import Alma
from glob import glob
import matplotlib.pyplot as plt
import xlsxwriter
import string
import os
from astropy.convolution import convolve, Gaussian2DKernel, Tophat2DKernel, Box2DKernel
from spectral_cube import SpectralCube
import warnings
warnings.filterwarnings('ignore')
from termcolor import colored
import pathlib
from math import floor


# Function to download observations
# Input: listobs,path,removecalibration, *specific
#   listobs (list<string,string>) : Each entry is an observations project code and sourcename.
#   path: (string) : Pathname of the directory of where to save files.
#   removecalibration(boolean): If true removes the calibration files.
#   specific*(string(s)): Optional argument to specify substring which the downloaded files should contain.
# Output: downloadedFiles
#   downloadedFiles (list<string>): List of the names of the downloaded files.
def downloadFiles(listObs, path, removecalibration = True, *specific):
    i = 0
    couldnotexpandproduct = []
    downloadedFiles = []
    for x in listObs:
        location = path + x[1]

        try:
            result = Alma.query_object('', project_code=x[0], source_name_alma = x[1])
        except:
            a=x[1].replace(" ","_")
            print(a)
            result = Alma.query_object('', project_code=x[0], source_name_alma=a)

        try:
            uids = np.unique(result['Member ous id'])
            link_list = Alma.stage_data(uids,expand_tarfiles=True)
            regex = re.compile('.*tar$')
            files = list(filter(regex.search, link_list['URL']))
            files = [file for file in files if x[1].upper() in file.upper()] #Remove Calibration
            if len(specific)>0:
                specific_files = []
                for sp in specific:
                    specific_files.append([file for file in files if sp.upper() in file.upper()])
                files = np.array(specific_files).flatten()

            if len(files) == 0:
                i += 1
                couldnotexpandproduct.append([x[0],x[1]])
            else:
                print(files)
                if not os.path.exists(location):
                    os.makedirs(location)
                Alma.cache_location = location
                try:
                    Alma.download_and_extract_files(files)
                    downloadedFiles.append(files)
                except Exception as err:
                    print(err)
        except Exception as err:
            print(err)
            print(x[1])
    downloadedFiles = np.array(downloadedFiles).flatten()
    return downloadedFiles


# Function to download observations
# Input: fileName,data
#   fileName: (string) : fileName of the excel document created.
#   data (table) : Table of observations from ASA
def excelTable(fileName, data):
    alphabet = string.ascii_uppercase[:]

    names = ('Project code', 'Source name', 'RA', 'Dec', 'Galactic longitude', 'Galactic latitude', 'Band',
             'Spatial resolution', 'Frequency resolution', 'Array', 'Mosaic', 'Integration', 'Release date',
             'Frequency support', 'Velocity resolution', 'Pol products', 'Observation date', 'PI name',
             'SB name', 'Proposal authors', 'Line sensitivity (10 km/s)', 'Continuum sensitivity', 'PWV',
             'Group ous id', 'Member ous id', 'Asdm uid', 'Project title', 'Project type', 'Scan intent',
             'Field of view', 'Largest angular scale', 'QA2 Status', 'COUNT', 'Science keyword',
             'Scientific category', 'ASA_PROJECT_CODE')

    workbook = xlsxwriter.Workbook(fileName)
    worksheet = workbook.add_worksheet()

    maxWidth = []
    for i in range(len(data[1])):
        maxWidth.append(0)

    for i in range(len(data) + 1):
        list = (names if i == 0 else data[i - 1])
        for j in range(len(list)):
            maxWidth[j] = max(min(max(maxWidth[j], len(str(list[j]))), 150), 6)

            pos = (str(alphabet[j]) + str(i + 1) if j < 26 else 'A' + str(alphabet[j - 26]) + str(i + 1))
            boldStatus = (True if i == 0 else False)
            wrapStatus = (True if i == 0 else False)
            bg_col = ('#F3F8EE' if i % 2 == 0 else '#DAECC9')
            format = workbook.add_format({'font': 'Times New Roman',
                                          'bold': boldStatus,
                                          'bg_color': bg_col,
                                          'border_color': 'white',
                                          'border': 1,
                                          'text_wrap': True,
                                          # 'shrink': True
                                          'valign': 'vcenter'
                                          })
            worksheet.write(pos, str(list[j]), format)
            worksheet.set_column(j, j, maxWidth[j] + 1)
    workbook.close()

# Function to query a table of observation from ASA
# Output: datalist,data
#   datalist(list<string,string>): List with each observations project code and source name
#   data(table): Table containin all queried observations

def initdatatable():
    data = Alma.query_object(''
                             , start_date='>01-01-2015'
                             , spatial_resolution='<0.1'
                             , integration_time='>1000'
                             , water_vapour='<2'
                             , band_list=['6'
            , '7'
            , '8'
            , '9'
                                          ]
                             , science_keyword=['Debris disks',
                                                'Disks around high-mass stars',
                                                'Disks around low-mass stars',
                                                'Astrochemistry',
                                                'HII regions',
                                                'High-mass star formation',
                                                'Low-mass star formation',
                                                'Inter-Stellar Medium (ISM)/Molecular clouds',
                                                'Intermediate-mass star formation',
                                                'Outflows, jets and ionized winds'
                                                ]
                             )

    cat = ['Active galaxies', 'Cosmology', 'Galaxy evolution']

    for x in cat:
        data = data[data['Scientific category'] != x]

    datalist = []

    for i in range(len(data)):
        datalist.append([data[i]['Project code'], data[i]['Source name']])

    return datalist,data

# Function to get the path for all fits file in a given directory
# Input: path
#   path(string): Path to search for fits files in.
# Output: result
#   result(list<string>): List with the paths for the files

def getAllFits(path):

    result = [y for x in os.walk(path) for y in glob(os.path.join(x[0], '*.fits'))]
    return result

def file_to_data(filename,info=False):
    from astropy.io import fits
    image_file = fits.open(filename)
    image_data = fits.getdata(filename)
    if info==True:
        print(image_data[0,:,:,:].shape)
        image_file.info()
    image_file.close()
    return image_data


def rms(data, return_rms=True, info=True):
    x_len = len(data[0])
    y_len = len(data[1])
    p = 0.03

    square_size = int(np.floor(np.multiply(x_len, p)))
    x = int(np.floor(np.divide(x_len, square_size)))
    y = int(np.floor(np.divide(y_len, square_size)))

    rms = np.inf  # inital value
    for i in range(0, x):
        for j in range(0, y):
            index_prel = [i, i + 1, j, j + 1]
            index_prel = np.multiply(square_size, index_prel)
            if np.isnan(data[index_prel[0]:index_prel[1], index_prel[2]:index_prel[3]]).any() or np.any(data[index_prel[0]:index_prel[1], index_prel[2]:index_prel[3]]==0):
                rms_prel = np.inf
            else:
                rms_prel = np.std(data[index_prel[0]:index_prel[1], index_prel[2]:index_prel[3]])
            if rms_prel < rms and rms_prel > 0:
                rms = rms_prel
                index = index_prel
    if info:
        print('Rms: ' + str(rms) + ' Jy/beam\nIndex for rms: ' + str(index))
    return (rms if return_rms else index)


def astro_plot(filename,zoom=-1,sigma=[3],color='purple',save=False):
    try:


        data_prel = file_to_data(filename)
        data = data_prel[0,0,:,:]
        zoom = astro_zoom(data,rms(data)*1000)

        fig_size = 4
        plt.figure(figsize=(fig_size*3,fig_size))
        plt.grid(True)
        plt.suptitle(filename.split('\\')[6],color='grey',size=14)

        plt.subplot(1,2,1)
        plt.imshow(data,cmap='hot')
        rectangle = plt.Rectangle((rms(data,False,False)[0],rms(data,False,False)[2])
                                  ,rms(data,False,False)[1]-rms(data,False,False)[0],
                                  rms(data,False,False)[3]-rms(data,False,False)[2],
                                  fill = None,
                                  ec = "white")
        plt.gca().add_patch(rectangle)
        plt.colorbar()
        plt.title(filename.split('\\')[6])
        ax = plt.axis()
        plt.axis((ax[0],ax[1],ax[3],ax[2]))

        plt.subplot(1,2,2)
        print('\n')
        sigma_plot(data,zoom,sigma,True,filename.split('\\')[6])


    except IndexError:
        print(colored('WARNING: IndexError: ' + filename,'blue'))
    except UnboundLocalError:
        print(colored('WARNING: IndexError: ' + filename,'blue'))


def sigma_plot(data, zoom, sigma=[3], subplot=True, title=''):
    data_full = data
    data_zoom = data_full[zoom[0]:zoom[1], zoom[2]:zoom[3]]
    max_value = np.max(data_zoom)  # maximum intensity
    contour_levels_max = max_value * 0.10 * range(1, 11)  # contour levels at 10%,20%,...,90% of max intensity
    std = 2*rms(data_full)  # standard deviation
    contour_levels_std_3 = sigma[0] * std * range(1, 5)  # contour levels at sigma

    if title == '':
        title = str(sigma[0]) + ' sigma'

    if not subplot:
        plt.figure()

    plt.grid(True)

    plt.title(title)
    plt.contour(data_zoom
                , levels=contour_levels_max
                , cmap='hot'
                )
    plt.colorbar()
    plt.contour(data_zoom
                , levels=contour_levels_std_3
                , colors='black'
                )

def astro_zoom(data,lim):
    data = (data[0, 0, :, :] if len(data.shape) == 4 else data[0, :, :] if len(data.shape) == 3 else data)
    data[np.where(np.isfinite(data[:, :]) == False)] = 0
    data_ravel_row = np.ravel(data)
    data_ravel_col = np.ravel(data, 'F')

    data_max = max(data_ravel_row)

    i = 0
    row_min = np.min(np.where(data_ravel_row > lim * data_max))
    row_max = np.max(np.where(data_ravel_row > lim * data_max))
    col_min = np.min(np.where(data_ravel_col > lim * data_max))
    col_max = np.max(np.where(data_ravel_col > lim * data_max))
    row = data.shape[0]
    col = data.shape[1]

    ind_row_low = floor(row_min / row)
    ind_col_low = row_min % row
    ind_row_high = floor(row_max / row)
    ind_col_high = row_max % row

    ind_row_low_2 = col_min % col
    ind_col_low_2 = floor(col_min / col)
    ind_row_high_2 = col_max % col
    ind_col_high_2 = floor(col_max / col)

    center = data.shape[0]/2 + data.shape[1]/2;

    ind = [np.zeros(4), np.zeros(4)]
    ind[0][0] = ind_row_low
    ind[1][0] = ind_col_low
    ind[0][1] = ind_row_high
    ind[1][1] = ind_col_high
    ind[0][2] = ind_row_low_2
    ind[1][2] = ind_col_low_2
    ind[0][3] = ind_row_high_2
    ind[1][3] = ind_col_high_2
    #return [int(floor((center - data.shape[0]/4))), int(floor((center + data.shape[0]/4))),int(floor((center - data.shape[0]/4))), int(floor((center + data.shape[0]/4)))]
    return [int(min(ind[:][0])), int(max(ind[:][0])), int(min(ind[:][1])), int(max(ind[:][1]))]