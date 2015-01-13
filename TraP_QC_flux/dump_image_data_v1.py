#!/usr/bin/env python
"""A specific example of how one might access the MonetDB database.

Try running with flag '--help' for information on the command line options.

"""

import optparse
import sys
import os
import tkp.config
import tkp.db
import csv

def dump_images(dbname,username,password,dataset_id, engine, host, port):
    tkp.db.Database(
        database=dbname, user=username, password=password,
        engine=engine, host=host, port=port
    )

    sources_query = """\
    SELECT  im.taustart_ts
            ,im.tau_time
            ,ex.f_int
            ,ex.f_int_err
            ,ax.xtrsrc
            ,rc.id as runcatid
            ,rc.dataset
            ,rc.wm_ra
            ,rc.wm_decl
            ,rf.avg_f_int
            ,rf.avg_f_int_sq
            ,im.freq_eff
            ,im.url
            ,sr.centre_ra
            ,sr.centre_decl
            ,im.rms_qc
            ,ex.extract_type
            ,im.rb_smaj
    FROM extractedsource ex
         ,assocxtrsource ax
         ,image im
         ,runningcatalog rc
         ,runningcatalog_flux rf
         ,skyregion sr
    WHERE rf.runcat = rc.id
      and ax.runcat = rc.id
      AND ax.xtrsrc = ex.id
      and ex.image = im.id
      AND rc.dataset = %s
      AND im.skyrgn=sr.id                                                                                                                                 
      ORDER BY rc.id
    """
    cursor = tkp.db.execute(sources_query, (dataset_id,))
    sources = tkp.db.generic.get_db_rows_as_dicts(cursor)
    print "Found", len(sources), "source datapoints"
    outfile_prefix = './ds_' + str(dataset_id) + '_'
    dump_list_of_dicts_to_csv(sources, outfile_prefix + 'sources_fluxQC.csv')

    return 0

def dump_list_of_dicts_to_csv(data, outfile):
    if data:
        with open(outfile, 'w') as fout:
            dw = csv.DictWriter(fout, fieldnames=sorted(data[0].keys()))
            dw.writerows(data)

