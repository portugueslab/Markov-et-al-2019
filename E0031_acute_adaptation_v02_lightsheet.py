import numpy as np
import pandas as pd

from stytra_config import ConfiguredStytra
from stytra.stimulation import Protocol
from stytra.stimulation.stimuli import CalibratingClosedLoop1D, Basic_CL_1D, GainChangerStimulus
from lightparam import Param
from stytra.stimulation.stimuli import GratingStimulus, AcuteClosedLoop1D


class ClosedLoop1DProtLS(Protocol):
    name = "E0031_acute_adaptation/v02_lightsheet"

    stytra_config = dict(
        tracking=dict(embedded=True, method="tail", estimator="vigor"),
        log_format="csv"
        )

    def __init__(self):
        super().__init__()

        self.grating_cycle = 10
        self.n_repeats_selfcalib = Param(10, limits=(0, 300))
        self.n_repeats_pre = Param(10, limits=(0, 300))
        self.n_repeats_exp = Param(40, limits=(0, 300))
        self.n_repeats_post = Param(10, limits=(0, 300))
        self.target_vel = Param(-18., limits=(-50, 50))
        self.calibrate_after = Param(10, limits=(0, 300))
        self.estim_gain = Param(-26., limits=(-100., 100.))

        inter_stim_pause = 7.5
        grating_vel = 10
        grating_duration = 15

        second_rev = 5
        first_rev = 2.5
        forward_duration = 0.35

        segments = [(0, first_rev),
                    (grating_vel, forward_duration),
                    (0, inter_stim_pause - first_rev - forward_duration),
                    (-grating_vel, grating_duration),
                    (0, second_rev),
                    (grating_vel, forward_duration),
                    (0, inter_stim_pause - second_rev - forward_duration)]

        t_base = [0]
        vel_base = [0]
        for s in segments:
            t_base.extend([t_base[-1], t_base[-1] + s[1]])
            vel_base.extend([s[0], ] * 2)

        self.t_base = np.array(t_base)
        self.vel_base = np.array(vel_base)

        self.ClosedLoop1DGratings = type("Stim", (Basic_CL_1D,
                                             GratingStimulus), {})

    def get_basic_block_df(self, n_reps):
        t = [0]
        vel = [0]
        for i in range(n_reps):
            t.extend(self.t_base + t[-1])
            vel.extend(self.vel_base)
        return pd.DataFrame(dict(t=t, base_vel=vel))

    def get_experimental_trials(self):
        df = self.get_basic_block_df(self.n_repeats_exp)

        ClosedLoop1DGratings = type("Stim", (AcuteClosedLoop1D,
                                             GratingStimulus), {})

        gains = [dict(change_to=dict(gain=g)) for g in [0, 1]] #0.66, 2]]
        #lags = [dict(change_to=dict(lag=l)) for l in [0.075, 0.225]]
        #drops_starts = [0., 0., 0.225, 0.075]
        #drops_ends = [0.075, 0.225, 0.3, 0.3]

        #drops = []
        #for start, end in zip(drops_starts, drops_ends):
        #    drops.append(dict(change_to=dict(gain_drop_start=start,
        #                                     gain_drop_end=end)))

        conditions = gains  # + lags + drops

        return ClosedLoop1DGratings(
            df_param=df,
            grating_angle=np.pi / 2,
            grating_period=self.grating_cycle,
            conditions_list=conditions)

    def get_stim_sequence(self):
        stimuli = []
        # # gratings

        # Calibration:
        stimuli.append(GainChangerStimulus(self.estim_gain))

        # pre trials:
        df = self.get_basic_block_df(self.n_repeats_pre)
        stimuli.append(
            self.ClosedLoop1DGratings(
                df_param=df,
                grating_angle=np.pi / 2,
                grating_period=self.grating_cycle,
                swimming_threshold=-2))

        stimuli.append(self.get_experimental_trials())

        return stimuli

if __name__ == "__main__":
    s = ConfiguredStytra(protocol=ClosedLoop1DProtLS())
