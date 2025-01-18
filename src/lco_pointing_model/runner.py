#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Author: José Sánchez-Gallego (gallegoj@uw.edu)
# @Date: 2024-07-31
# @Filename: runner.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from __future__ import annotations

import asyncio
import pathlib

import polars
from clu.legacy import TronConnection
from clu.legacy.types.parser import Reply
from pydantic import BaseModel, Field, create_model

from sdsstools import get_sjd

from lco_pointing_model import log
from lco_pointing_model.sample import get_random_sample, to_icrs


SCHEMA: dict[str, polars.DataTypeClass] = {
    "alt": polars.Float64,
    "az": polars.Float64,
    "rot": polars.Float64,
    "ra": polars.Float64,
    "dec": polars.Float64,
    "done": polars.Boolean,
    "failed": polars.Boolean,
    "mjd": polars.Int32,
    "gimg": polars.Int32,
    "n_cameras_solved": polars.Int32,
    "ra_bore": polars.Float64,
    "dec_bore": polars.Float64,
    "offset_ra": polars.Float64,
    "offset_dec": polars.Float64,
    "offset_rot": polars.Float64,
    "tai_ref": polars.Float64,
    "ra_tel": polars.String,
    "dec_tel": polars.String,
    "st_tel": polars.String
}


class PointingDataBase(BaseModel):
    """Pointing model data."""

    alt: float
    az: float
    rot: float
    ra: float | None = None
    dec: float | None = None
    done: bool = False
    failed: bool = False
    mjd: float = Field(default=get_sjd("LCO"))
    gimg: int | None = None
    n_cameras_solved: int = 0
    tai_ref: float | None = None
    ra_bore: float | None = None
    dec_bore: float | None = None
    offset_ra: float | None = None
    offset_dec: float | None = None
    offset_rot: float | None = None
    ra_tel: str | None = None
    dec_tel: str | None = None
    st_tel: str | None = None


async def query_tcs(waittime=5):
    await asyncio.sleep(waittime)

    reader, writer = await asyncio.open_connection(
        'c100tcs.lco.cl', 4242)

    raCmd = "alp\r\n"
    writer.write(raCmd.encode())
    await writer.drain()

    data = await reader.read(100)
    ra = data.decode().strip()

    decCmd = "del\r\n"
    writer.write(decCmd.encode())
    await writer.drain()

    data = await reader.read(100)
    dec = data.decode().strip()

    stCmd = "st\r\n"
    writer.write(stCmd.encode())
    await writer.drain()

    data = await reader.read(100)
    st = data.decode().strip()

    print(ra)
    print(dec)
    print(st)

    writer.close()
    await writer.wait_closed()
    return ra, dec, st


async def get_pointing_data(
    npoints: int | None,
    output_file: str | pathlib.Path,
    reuse_file: bool = True,
    overwrite: bool = False,
    write_csv: bool = True,
    alt_range: tuple[float, float] = (30, 90),
    az_range: tuple[float, float] = (0, 359),
    write_log: bool = True,
):
    """Collects pointing model data.

    Parameters
    ----------
    npoints
        The number of grid points to collect.
    output_file
        The output file. If the file exists and ``overwrite=False``, the grid will
        be read from the file. The file will be written in Parquet format so it is
        expected to have ``.parquet`` extension.
    reuse_file
        If True, will reuse the file if it exists and will not generate a new grid.
        If False and the file exists, it will raise an error unless ``overwrite=True``.
    overwrite
        If True, overwrites the output file if it exists.
    write_csv
        If True, writes the output to a CSV file as well.
    alt_range
        The range of altitudes to sample.
    az_range
        The range of azimuths to sample.
    write_log
        If True, writes a log file with the same name as the output file but with
        the extension replaced to ``.log``

    Returns
    -------
    df
        A Polars DataFrame with the pointing model data.

    """

    PointingData = create_model(
        "PointingData",
        __base__=PointingDataBase,
    )

    output_file = pathlib.Path(output_file)

    if write_log:
        log.start_file_logger(str(output_file.with_suffix(".log")), rotating=False)

    ### Recover pointing grid or create a new one. ###

    if output_file.exists() and reuse_file:
        log.warning(f"Found file {output_file!s}. NOT generating a new grid.")
        data_df = polars.read_parquet(output_file)
        data = [PointingData(**row) for row in data_df.to_dicts()]
    else:
        if output_file.exists() and overwrite is False:
            raise FileExistsError("output_file exists and overwrite=False.")

        if npoints is None:
            raise ValueError("npoints is required if output_file does not exist.")

        log.info("Generating new pointing model grid.")
        altaz = get_random_sample(npoints, alt_range=alt_range, az_range=az_range)
        data: list[PointingDataBase] = []
        for alt, az in altaz:
            data.append(PointingData(alt=alt, az=az, rot=0.0))

        write_data(data, output_file, write_csv=write_csv)

    ### Create Tron connection and instantiate TCC helper. ###

    log.info("Creating connection to Tron and waiting for keys.")
    tron = TronConnection(
        "LCO.pointing_model",
        "sdss5-hub.lco.cl",
        models=["lcotcc"],
    )
    await tron.start()
    await asyncio.sleep(5)

    ### Loop over the data and solve the fields. ###

    for irow, pdata in enumerate(data):
        if pdata.done or pdata.failed:
            log.warning(f"Skipping row {irow!r} as it is already done or failed.")
            continue

        alt = pdata.alt
        az = pdata.az

        ra, dec = to_icrs(alt, az)
        pdata.ra = ra
        pdata.dec = dec

        ### Slew to the target. ###

        log.warning(f"Slewing to target #{irow + 1}: RA={ra:.4f}, Dec={dec:.4f}.")

        slew_cmd = await tron.send_command(
            "lcotcc",
            f"target {ra}, {dec} /posAngle=270.0",
            callback=log_reply,
            time_limit=600,
        )
        if slew_cmd.status.did_fail:
            log.error("Failed to slew to target. Skipping this pointing.")
            pdata.failed = True
            continue

        # Make sure the guider knows the position of the field.
        log.debug("Setting a fake jaeger field.")
        await tron.send_command(
            "jaeger",
            f"configuration fake-field {ra} {dec} 270.0",
        )

        ### Solve the field. ###

        log.info("Astrometrically solving the field with cherno.")
        exp_time: float = 5

        while True:
            # cmd_time = await tron.send_command("tcc", "show time")
            # tai0 = cmd_time.replies.get("tai")[0]
            tai0 = 5
            print("hacking tai0", 5)

            L = await asyncio.gather(
                query_tcs(exp_time/2.),
                tron.send_command(
                    "cherno",
                    f"acquire -t {exp_time} --mode astrometrynet --no-continuous --no-apply",
                    callback=log_reply,
                )

            )

            tcs_data, cmd_acq = L

            print("tcs_data", tcs_data)

            did_fail = cmd_acq.status.did_fail
            acquisition_valid = cmd_acq.replies.get("acquisition_valid")[0]

            try:
                astrometry_fit = cmd_acq.replies.get("astrometry_fit")
                if float(astrometry_fit[2]) == -999.0:
                    log.warning(
                        "Field was not solved with fit_SP or failed to acquire."
                    )
                    astrometry_fit = None
            except KeyError:
                log.warning("The astrometry_fit keyword was not output.")
                astrometry_fit = None

            if did_fail or astrometry_fit is None or not acquisition_valid:
                exp_time += 10
                if exp_time <= 25:
                    log.error(
                        "Failed to solve field. Retrying "
                        f"with exp_time={exp_time} seconds."
                    )
                    continue

                log.error("Failed to solve field. Skipping this pointing.")
                pdata.failed = True
                break

            pdata.gimg = int(astrometry_fit[0])
            pdata.n_cameras_solved = int(astrometry_fit[1])
            pdata.ra_bore = float(astrometry_fit[2])
            pdata.dec_bore = float(astrometry_fit[3])
            pdata.offset_ra = float(astrometry_fit[7])
            pdata.offset_dec = float(astrometry_fit[8])
            pdata.offset_rot = float(astrometry_fit[9])
            pdata.tai_ref = float(tai0 + exp_time / 2) # not used?
            pdata.ra_tel = str(tcs_data[0])
            pdata.dec_tel = str(tcs_data[1])
            pdata.st_tel = str(tcs_data[2])

            # Override MJD in case the pointing is done at some other point in time
            pdata.mjd = get_sjd("LCO")

            ### Run ptcorr to get the pointing corrections. ###

            log.info("Retrieving pdata.")
            pdata.done = True
            break
            # TODO

        ### Update data and write to disk. ###
        write_data(data, output_file, write_csv=write_csv)

    return data


def write_data(
    data: list[PointingDataBase],
    output_file: pathlib.Path,
    write_csv: bool = True,
):
    """Writes the data to disk.

    Parameters
    ----------
    data
        The list of pointing data.
    output_file
        The output file. Expected to be the path to a Parquet file.
    write_csv
        If True, writes the data to a CSV file as well. The file will have the same
        name as ``output_file`` but with the extension replaced to ``.csv``.

    """

    schema = SCHEMA.copy()

    data_df = polars.DataFrame([d.model_dump() for d in data], schema=schema)
    data_df.write_parquet(output_file)

    if write_csv:
        csv_file = output_file.with_suffix(".csv")
        data_df.write_csv(csv_file)


def log_reply(reply: Reply):
    """Logs a reply to a command."""

    log.debug(reply.string)
