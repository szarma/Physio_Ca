import pandas as pd
import re
import numpy as np
import warnings
num_alph = re.compile("\d[a-zA-Z]")

class Protocol(pd.DataFrame):
    def add_conc_and_unit_columns(self):
        cols = list(self.columns)
        for i in self.index:
            concentration = self.loc[i, "concentration"].replace(" ", "")
            self.loc[i, "concentration"] = concentration
            matches = num_alph.findall(concentration)
            if len(matches) > 1:
                raise ImportError("at least one concentration does not have the appropriate form (quantity, unit)")
            elif len(matches) == 0:
                if concentration[0].isdigit():
                    conc = concentration
                    unit = ""
                else:
                    conc = ""
                    unit = concentration

            else:
                isplit = num_alph.search(concentration).start() + 1
                conc = concentration[:isplit]
                unit = concentration[isplit:]
            self.loc[i, "conc"] = conc
            self.loc[i, "unit"] = unit
        self = self[cols[:2] + ["conc", "unit"] + cols[2:]]

    @classmethod
    def from_df(cls, df):
        protocol = cls(df)
        for c in ["conc","unit"]:
            if c in df.columns:
                del df[c]
        protocol.add_conc_and_unit_columns()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category = UserWarning)
            protocol.legs = protocol.get_legs()
        return protocol

    @classmethod
    def from_file(cls, path, tend, t0=0):
        lines = open(path).readlines()
        for il,line in enumerate(lines):
            n_commas = line.count(",")
            if n_commas!=3:
                raise ValueError(f"Line {il} of {path} has {n_commas}, but each line must have exactly 3.")
        df = pd.read_csv(path, na_values = [" "*j for j in range(1,10)])
        if len(df)==0:
            raise pd.errors.EmptyDataError
        beginNa = df.index[df["begin"].isna()]
        df.loc[beginNa, "begin"] = "00:%02i" % t0
        df["t_begin"] = [pd.Timedelta("00:" + v).total_seconds() for v in df["begin"]]
        df["t_end"] = [pd.Timedelta("00:" + v).total_seconds() if isinstance(v, str) else np.nan for v in
                             df["end"]]
        if df["t_end"].isna().any():
            endNa = df.index[df["t_end"].isna()]
            df.loc[endNa, "t_end"] = max(tend, df["t_end"].max() )
        protocol = cls.from_df(df)
        protocol.path_to_file = path
        return protocol

    def get_legs(self, parse_output=False, verbose=False):
        out = {}
        timepoints = set(self[["t_begin", "t_end"]].values.flatten())
        timepoints = sorted(timepoints)
        for t in timepoints:
            dft = self.query(f"t_begin<={t} and {t}<t_end").copy()
            if verbose:
                print (t, "\n", dft,"\n","="*5)

            if len(dft["compound"].unique()) != len(dft):
                raise ValueError("Inconsistent treatments. Legs that should not overlap, do overlap in time. Edit the file and try again.")
            if parse_output:
                v = {row["compound"].lower(): (row["conc"], row["unit"]) for ir, row in dft.iterrows()}
            else:
                v = {row["compound"].lower(): row["concentration"] for ir, row in dft.iterrows()}
            if len(v):
                out[t, dft["t_end"].min()] = v
            for k in out:
                if "glucose" not in out[k]:
                    out[k]["glucose"] = ("0","")
        return out

    def get_scheme(self, symbolDict):
        prtl = []
        legs = self.get_legs(parse_output = True)

        for k in legs:
            leg = legs[k]
            letter = leg["glucose"][0]
            for comp in leg:
                if comp == "glucose":
                    continue
                symbol_key = "%s[%s%s]" % (comp, leg[comp][0], str(leg[comp][1]))
                if symbol_key not in symbolDict:
                    newk = "*" * len(symbolDict) if len(symbolDict) < 3 else "¹²³⁴⁵⁶⁷⁸⁹⁰"[len(symbolDict) - 3]
                    symbolDict[symbol_key] = newk
                letter += symbolDict[symbol_key]
            prtl += [letter]
        return "-".join(prtl)
    
    def copy(self, deep=True):
        protocol = super(Protocol, self).copy(deep=deep)
        protocol = self.from_df(protocol)
        return protocol
        
    def plot_protocol( self, ax, color=(.95,)*3, hspace = 0.2, only_number = True, linenkwargs={}, label = "close", fontsize=10, offset=.05):
        if only_number:
            col = "conc"
        else:
            col = "concentration"

        fig = ax.get_figure()
        compounds = list(self["compound"].unique())
        # put glucose first
        for ic,comp in enumerate(compounds):
            if "glc" in comp.lower() or "glu" in comp.lower() or comp=="" :
                compounds.pop(ic)
                compounds = [comp] + compounds
                break
        nComp = len(compounds)
        figheight = fig.get_figheight()
        figheight = figheight*72
        axHeight = nComp*fontsize*(1+2*hspace)/figheight
        pos = ax.get_position()
        axp = fig.add_axes([pos.x0, pos.y0+pos.height*(1+offset), pos.width, axHeight], label="protocol")
        axp.set_ylim(0, nComp)
        fig.canvas.draw()
        for ic,comp in enumerate(compounds):
            df = self.query(f"compound=='{comp}'")
            for ii in df.index:
                t0,t1 = df.loc[ii,["t_begin","t_end"]]
                conc = df.loc[ii,col]
                x = [t0,t1,t1,t0,t0]
                y = np.array([1,1,0,0,1])*(1-hspace)+ic+hspace/2
                c = df.loc[ii].get("color",color)
                if "edgecolor" not in linenkwargs:
                    linenkwargs["edgecolor"] = "black"

                if "linewidth" not in linenkwargs:
                    linenkwargs["linewidth"] = 1

                ln = axp.fill(x,y,color=c,**linenkwargs)[0]
                ln.set_clip_on(False)
                ftc = df.loc[ii].get("fontcolor","black")
                yt = ic+.5 - fontsize/figheight#y[:-1].mean()
                axp.text(t0,yt, " "+str(conc),va="center", ha="left", color=ftc, size = fontsize)
                #ax.plot(x,y,color="black",lw=1)
            if label == "close":
                axp.text(df.t_begin.min(),yt,comp+" ",va="center", ha="right", fontsize = fontsize)
            if label == "left":
                axp.text(self.t_begin.min(),yt,comp+" ",va="center", ha="right", fontsize = fontsize)
        axp.set_clip_on(False)
        #ax.set_yticks([])
        #mystyle_axes(ax)
        return axp
