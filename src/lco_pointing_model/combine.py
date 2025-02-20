#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: Conor Sayres (csayres@uw.edu)
# @Date: 2024-08-29
# @Filename: combine.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import datetime
from astropy.coordinates import SkyCoord

import polars


def getTpointHeader(utcDateStr):
    """Returns list of lines to include in TPOINT header"""
    currDateStr = datetime.datetime.now().isoformat()
    # header = [
    #     "! Caption record:",
    #     "SDSS APO 2.5m Pointing Data " + currDateStr,
    #     "! Option record: ALTAZ followed by rotator code, if applicable",
    #     ": ALTAZ",
    #     "! Run parameters: telescope latitude deg min sec",
    #     "32 46 49.3",
    #     "! Pointing data (TPOINT format #4):",
    #     "!",
    #     "! All angles are in degrees",
    #     "! Azimuth uses the convention: N=0, E=90 (unlike the TCC, which uses S=0, E=90)",
    #     "!",
    #     "!  Desired Phys            Actual Mount          Rot Phys",
    #     "! Az          Alt         Az          Alt",
    # ]
    header = [
        "duPont " + currDateStr + " SDSS-V FPS",
        ": J2000",
        ": EQUAT",
        ": NODA",
        ": ALLSKY",
        "! Latitude, UTC date, temp, pressure, height in m, relative humidity",
        "-29 00 26.8 %s 13  777  2282.0 0.06"%utcDateStr,
        "",
        "!RAtarg Dectarg ppm1 ppm2 epoch RA DEC sth stm",
    ]
    return header


def processFiles(fileInputList, fileOutput, tpointOn=False):
    """Get a TPOINT formatted file from one or more pointing data csv files

    Inputs:
    - file_input_list: list of paths to csv files to be used
    - file_output: path to output file

    Returns None
    """
    dfList = []
    for fin in fileInputList:
        dfList.append(polars.read_csv(fin))
    df = polars.concat(dfList)
    nMeas = len(df)
    df = df.filter(~polars.col("failed"))
    df = df.filter(polars.col("done"))
    nGood = len(df)
    print("using %i of %i pointings" % (nGood, nMeas))

    sky_coords = SkyCoord(df["ra_bore"].to_numpy(), df["dec_bore"].to_numpy(), frame="icrs", unit="deg")
    sky_coords = sky_coords.to_string("hmsdms", sep=" ")
    ppm1 = [0]*len(sky_coords)
    ppm2 = [0]*len(sky_coords)
    epoch = [2000.0]*len(sky_coords)
    if tpointOn:
        RA_tel = df["ra_tel_tp"]
        Dec_tel = df["dec_tel_tp"]
    else:
        RA_tel = df["ra_tel"]
        Dec_tel = df["dec_tel"]

    RA_tel = [x.replace(":", " ") for x in RA_tel]
    Dec_tel = [x.replace(":", " ") for x in Dec_tel]
    ST_tel = []
    for st in df["st_tel"]:
        h,m,s = st.split(":")
        dm = float(m) + float(s)/60.
        ST_tel.append("%s %.2f"%(h, dm))

    # convert from tcc az convention to TPOINT az convention (N=0, E=90)

    fileLines = getTpointHeader(df["ut_date"][0])



    for ii in range(nGood):
        fileLines.append(
            "%s %.1f %.1f %.1f %s %s %s"
            % (
                sky_coords[ii].replace("+", ""), # dont write+ for +dec
                ppm1[ii],
                ppm2[ii],
                epoch[ii],
                RA_tel[ii],
                Dec_tel[ii],
                ST_tel[ii]
            )
        )

    print("writing", fileOutput)
    with open(fileOutput, "w") as f:
        for line in fileLines:
            f.write(line + "\n")

# if __name__ == "__main__":
#     processFiles(["../../mjd60723.csv"], "../../mjd60723_tpoint_on.dat")
