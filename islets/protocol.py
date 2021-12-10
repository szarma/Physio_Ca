import pandas as pd
import re
import numpy as np
import warnings


class Protocol(pd.DataFrame):
    @classmethod
    def from_file(cls, path, tend, t0=0):
        num_alph = re.compile("\d\D")
        protocol = pd.read_csv(path)
        protocol["begin"][protocol["begin"].isna()] = "00:%02i" % t0
        protocol["t_begin"] = [pd.Timedelta("00:" + v).total_seconds() for v in protocol["begin"]]
        protocol["t_end"] = [pd.Timedelta("00:" + v).total_seconds() if isinstance(v, str) else np.nan for v in
                             protocol["end"]]
        if protocol["t_end"].isna().any():
            protocol.loc[protocol["t_end"].isna().idxmax(), "t_end"] = max(tend, protocol["t_end"].max() + 1)
        cols = list(protocol.columns)
        for i in protocol.index:
            concentration = protocol.loc[i, "concentration"].replace(" ", "")
            protocol.loc[i, "concentration"] = concentration
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
            protocol.loc[i, "conc"] = conc
            protocol.loc[i, "unit"] = unit
        protocol = protocol[cols[:2] + ["conc", "unit"] + cols[2:]]
        protocol = cls(protocol)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            protocol.legs = protocol.get_legs()
        return protocol

    def get_legs(self, parse_output=False):
        out = {}
        for t in self["t_begin"].sort_values():
            dft = self.query(f"t_begin<={t} and {t}<t_end").copy()
            if len(dft["compound"].unique()) != len(dft):
                raise ValueError("Inconsistent treatments.")
            if parse_output:
                v = {row["compound"]: (row["conc"], row["unit"]) for ir, row in dft.iterrows()}
            else:
                v = {row["compound"]: row["concentration"] for ir, row in dft.iterrows()}
            out[t, dft["t_end"].min()] = v
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
