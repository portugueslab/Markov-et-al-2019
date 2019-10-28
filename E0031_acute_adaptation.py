import numpy as np
import pandas as pd

from stytra_config import ConfiguredStytra
from stytra_config.protocols.omr.E0030_long_term_adaptation_v01_normal import ClosedLoop1DProt
from stytra.stimulation.stimuli import GratingStimulus, AcuteClosedLoop1D


class AcuteClosedLoop1DProt(ClosedLoop1DProt):
    name = "E0031_acute_adaptation"

    def get_experimental_trials(self):
        t = [0]
        vel = [0]
        for i in range(self.n_repeats_exp):
            t.extend(self.t_base + t[-1])
            vel.extend(self.vel_base)
        df = pd.DataFrame(dict(t=t, base_vel=vel))

        ClosedLoop1DGratings = type("Stim", (AcuteClosedLoop1D,
                                             GratingStimulus), {})

        gains = [dict(change_to=dict(gain=g)) for g in [0, 1, 0.33, 0.66, 1.33, 1.66, 2]]
        lags = [dict(change_to=dict(lag=l)) for l in [0.075, 0.15, 0.225, 0.3]]
        # shunted_lag = [dict(change_to=dict(lag=0.5, shunted=True))]
        drops_starts = [0., 0., 0., 0., 0.225, 0.15, 0.075]
        drops_ends = [0.075, 0.15, 0.225, 0.3, 0.3, 0.3, 0.3]

        drops = []
        for start, end in zip(drops_starts, drops_ends):
            drops.append(dict(change_to=dict(gain_drop_start=start,
                                                  gain_drop_end=end)))

        conditions = gains + lags + drops

        return ClosedLoop1DGratings(
                df_param=df,
                grating_angle=np.pi / 2,
                grating_period=self.grating_cycle,
                conditions_list=conditions)


if __name__ == "__main__":
    s = ConfiguredStytra(protocol=AcuteClosedLoop1DProt())
