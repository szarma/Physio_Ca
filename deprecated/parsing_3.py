import pandas as pd

saveFilename = "/data/.tmp/protocols.csv"
    
def update_protocols(sheet=2):
    if sheet==1:
        sandras_protocols_url = "https://docs.google.com/spreadsheets/d/1AkF1NzddXQ_ATdCvvqX_eXET6kHRMwcoExxHZBY9Quw/export?gid=0&format=csv"
    if sheet==2:
        sandras_protocols_url = "https://docs.google.com/spreadsheets/d/1AkF1NzddXQ_ATdCvvqX_eXET6kHRMwcoExxHZBY9Quw/export?format=csv&id=1AkF1NzddXQ_ATdCvvqX_eXET6kHRMwcoExxHZBY9Quw&gid=1663512963"
    from urllib.request import urlretrieve
    from os.path import expanduser
    urlretrieve(sandras_protocols_url,saveFilename)
    return None

def import_protocols(parseMetadata=False,user="Sandra", verbose=False):
    from os.path import expanduser
    protocols = pd.read_csv(saveFilename)
    for col in protocols:
        if protocols[col].isna().all():
            del protocols[col]
    for col in protocols:
        if col[-2]=="." and col[:-2] not in protocols:
            protocols[col[:-2]] = protocols[col]
            del protocols[col]
    protocols = protocols.sort_values("date")
    protocols.index = range(len(protocols))
    protocols.loc[protocols.query("microscope=='CF_SP8'").index,"microscope"] = 'leica CF_SP8'
    for col in protocols:
        try:
            protocols[col] = protocols[col].apply(str.strip)
        except:
            pass
    if not parseMetadata:
        return protocols
    
    from islets.Recording import Recording
    import os 
    collect = []
    for cd, dirs, fs in os.walk("/data/"+user):
        for fn in fs:
            if fn.endswith(".lif") or fn.endswith(".nd2"):
                arow = {
                    "path": os.path.join(cd,fn),
                    "filename": fn[:-4].lower(),
                    "ext": fn[-3:]
                       }
                collect += [arow]

    localExps = pd.DataFrame(collect)
    df = pd.merge(protocols,localExps,how="left",left_on="experiment",right_on="filename")
    attribs = ['SizeT', 'SizeX', 'SizeY', 'pxSize', 'bit depth',
               'Frequency', 'Start time', 'End time', 'Duration',
#                'frame_range'
              ]
    for attr in attribs: df[attr] = None
    df["Parsing Error"] = [""]*len(df)
    for exp, df_exp in df.groupby("experiment"):
        path = df_exp.path.iloc[0]
        assert (df_exp.path == path).all()
        rec = Recording(path)
        rec.calc_gaps()
        for ix,row in df_exp.iterrows():
            try:
                rec.import_series(row.series, onlyMeta=True)
            except:
#                 from sys import exc_info
                if verbose:
                    print (f"could not import metadata for {row.series} for {row.experiment}. It contains the following series with the following frequencies: ")
                    print (rec.metadata[["Name","Frequency"]])
                df.loc[ix,"Parsing Error"] = f"Error: Recording contains the following series:\n"+\
                    str(rec.metadata[["Name","Frequency",'SizeT', 'SizeX', 'SizeY', 'pxSize',"gap"]])
#                 print (f"could not import metadata for {row.series} from {path}.")
#                 print ("Error:")
#                 print (exc_info())
                continue
            tmpMeta = rec.Series[row.series]["metadata"]
            df.loc[ix,attribs] = tmpMeta[attribs]
    return df