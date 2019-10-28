import numpy as np
import pandas as pd

from stytra_config import ConfiguredStytra
from stytra_config.protocols.omr.E0030_long_term_adaptation_v01_normal import ClosedLoop1DProt
from stytra.stimulation.stimuli import GratingStimulus, GainLagClosedLoop1D


class OpenLoopClosedLoop1DProt(ClosedLoop1DProt):
    name = "E0030_long_term_adaptation/v02_lag"

    def get_experimental_trials(self):
        t = [0]
        vel = [0]
        for i in range(self.n_repeats_exp):
            t.extend(self.t_base + t[-1])
            vel.extend(self.vel_base)
        df = pd.DataFrame(dict(t=t, base_vel=vel))

        ClosedLoop1DGratings = type("Stim", (GainLagClosedLoop1D,
                                             GratingStimulus), {})

        return ClosedLoop1DGratings(
                df_param=df,
                grating_angle=np.pi / 2,
                grating_period=self.grating_cycle,
                lag=0.225)


if __name__ == "__main__":
    s = ConfiguredStytra(protocol=OpenLoopClosedLoop1DProt())