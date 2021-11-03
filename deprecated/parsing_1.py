import pandas as pd

# def readSingleMetadata(series,seriesDir):
#     from os import listdir
#     from os.path import sep
#     fs = [f for f in listdir(seriesDir) if series in f and "txt" in f]
#     if len(fs)==0:
#         print (f"could not find {series} metadata file in {seriesDir}.\n\tCheck manually if you need, I'll continue.")
#         return False
#     if len(fs)>1:
#         print(f" found two files that have the form indicating they contain metadata for {f}. I'll take the first one, but beware...")
#     f = fs[0]
#     add = pd.read_csv(seriesDir+f).iloc[0]
# #     add["filename"] = f.split(".")[0]
#     add["filename"] = sep.join(seriesDir[:-1].split(sep)[3:]+[f])
#     return add


def update_protocols(sheet=1):
    if sheet==1:
        sandras_protocols_url = "https://docs.google.com/spreadsheets/d/1AkF1NzddXQ_ATdCvvqX_eXET6kHRMwcoExxHZBY9Quw/export?gid=0&format=csv"
    if sheet==2:
        sandras_protocols_url = "https://docs.google.com/spreadsheets/d/1AkF1NzddXQ_ATdCvvqX_eXET6kHRMwcoExxHZBY9Quw/export?format=csv&id=1AkF1NzddXQ_ATdCvvqX_eXET6kHRMwcoExxHZBY9Quw&gid=1663512963"
    from urllib.request import urlretrieve
    from os.path import expanduser
    urlretrieve(sandras_protocols_url,expanduser("~/protocols.csv"))
    return None

def import_protocols(
    parseMetadata=True,
    colToNum=["pH","glucose"],
    colToDel=[
        "mouse","sex",   # because not useful for now
        "frequency","resolution" # because can be read from metadata and is more reliable
    ]):
    from os.path import expanduser
    protocols = pd.read_csv(expanduser("~/protocols.csv"))
    for col in colToNum:
        protocols[col+"_"]  = protocols[col].copy()
        protocols[col] = pd.to_numeric(protocols[col],errors="coerce")
    protocols.treatment = protocols.treatment.astype(str)
    for col in colToDel:
        try:
            del protocols[col]
            print ("deleted ", col, "as not useful")
        except:
            pass
    protocols = protocols.sort_values("date")
    print ("when pH is NAN I will assume it is neutral (7.4).")
    protocols.loc[protocols.pH.isna(),"pH"] = 7.4
    protocols.index = range(len(protocols))
    protocols = protocols[[c for c in protocols if c!="magnification"]+["magnification"]]
    protocols.loc[protocols.query("microscope=='CF_SP8'").index,"microscope"] = 'leica CF_SP8'
    if not parseMetadata: return protocols
    from os.path import isdir
    mentioned=[]
    for ix,row in protocols.iterrows():
        isNikon = "leica" not in row.microscope.lower()
        if isNikon:
            ser = row.experiment.split(".")[0]
            serDir = expanduser("~/local_data")+f"/Sandra/NIKON/{row.date}/{ser}/"
        else:
            exp, ser = row.experiment.split("_")
            ser = ser.split()[0]
            serDir = expanduser("~/local_data")+f"/Sandra/{row.date}/{exp}/"
            
        if isdir(serDir):
            if "-" not in ser:
                add = readSingleMetadata(ser,serDir)
                if add is False:
                    continue
            else:
                if isNikon:
                    print ("stitching not implemented for ", row.experiment, "... skipping.")
                    continue
                serSet = [int(el) for el in ser.replace("Series","").split("-")]
                serSet[1] += 1
                serSet = ["Series%03i"%i for i in range(*serSet)]
                adds = [readSingleMetadata(ser,serDir) for ser in serSet]
                adds = [add for add in adds if add is not False]
                if len(adds)==0:
                    continue
                adds = pd.DataFrame(adds)
                Ymode = int(adds.Y.mode())
                if Ymode==8000: ## line scans
                    adds = adds.query(f"Y=={Ymode}")
                    frmed = adds.freq.median()
                    adds = adds.query(f"freq<{frmed*1.01} and freq>{frmed/1.01}")
                assert adds.X.std()==0
                assert adds.Y.std()==0
                assert adds.freq.std()/adds.freq.mean()<1e-3
                add = adds.iloc[0].copy()
                add["T"] = adds["T"].sum()
                add.freq = adds.freq.mean()
                add.filename = "\n".join(adds.filename)

            for k in add.keys():
                protocols.loc[ix,k] = add[k]
        else:
            if serDir not in mentioned:
                mentioned += [serDir]
                print(f"{serDir} does not exist.")
    try: del protocols["Unnamed: 0"]
    except: pass
    protocols["Tmax"] = protocols["T"]/protocols.freq
    return protocols