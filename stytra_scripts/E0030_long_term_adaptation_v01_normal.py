import numpy as np
import pandas as pd

from stytra import Stytra
from stytra.stimulation import Protocol
from stytra.stimulation.stimuli import CalibratingClosedLoop1D, GratingStimulus, Basic_CL_1D
from lightparam import Param


class ClosedLoop1DProt(Protocol):
    name = "E0030_long_term_adaptation/v01_normal"

    stytra_config = dict(
        tracking=dict(embedded=True, method="tail", estimator="vigor"),
        log_format="csv"
        )

    def __init__(self):
        super().__init__()

        self.grating_cycle = 10
        self.n_repeats_selfcalib = Param(10, limits=(0, 300))
        self.n_repeats_pre = Param(10, limits=(0, 300))
        self.n_repeats_exp = Param(210, limits=(0, 300))
        self.n_repeats_post = Param(10, limits=(0, 300))
        self.target_vel = Param(-18., limits=(-50, 50))
        self.calibrate_after = Param(10, limits=(0, 300))

        inter_stim_pause = 7.5
        grating_vel = 10
        grating_duration = 15

        self.t_base = np.array([0, inter_stim_pause,
                                inter_stim_pause, inter_stim_pause + grating_duration,
                                inter_stim_pause + grating_duration, 2 * inter_stim_pause +grating_duration])
        self.vel_base = np.array([0, 0, -grating_vel, -grating_vel, 0, 0])

        self.ClosedLoop1DGratings = type("Stim", (Basic_CL_1D,
                                             GratingStimulus), {})

    def get_experimental_trials(self):
        t = [0]
        vel = [0]
        for i in range(self.n_repeats_exp):
            t.extend(self.t_base + t[-1])
            vel.extend(self.vel_base)
        df = pd.DataFrame(dict(t=t, base_vel=vel))

        return self.ClosedLoop1DGratings(
                df_param=df,
                grating_angle=np.pi / 2,
                grating_period=self.grating_cycle,
                swimming_threshold=-2)

    def get_stim_sequence(self):
        stimuli = []
        # # gratings

        t = [0]
        vel = [0]
        for i in range(self.n_repeats_selfcalib):
            t.extend(self.t_base + t[-1])
            vel.extend(self.vel_base)
        df = pd.DataFrame(dict(t=t, base_vel=vel))

        ClosedLoop1DCalibratingGratings = type("Stim", (CalibratingClosedLoop1D,
                                             GratingStimulus), {})

        stimuli.append(
            ClosedLoop1DCalibratingGratings(
                df_param=df,
                grating_angle=np.pi / 2,
                grating_period=self.grating_cycle,
                swimming_threshold=-2,
                target_avg_fish_vel=self.target_vel,
                calibrate_after=self.calibrate_after))

        t = [0]
        vel = [0]
        for i in range(self.n_repeats_pre):
            t.extend(self.t_base + t[-1])
            vel.extend(self.vel_base)
        df = pd.DataFrame(dict(t=t, base_vel=vel))

        stimuli.append(
            self.ClosedLoop1DGratings(
                df_param=df,
                grating_angle=np.pi / 2,
                grating_period=self.grating_cycle,
                swimming_threshold=-2))

        stimuli.append(self.get_experimental_trials())

        t = [0]
        vel = [0]
        for i in range(self.n_repeats_post):
            t.extend(self.t_base + t[-1])
            vel.extend(self.vel_base)
        df = pd.DataFrame(dict(t=t, base_vel=vel))

        stimuli.append(
            self.ClosedLoop1DGratings(
                df_param=df,
                grating_angle=np.pi / 2,
                grating_period=self.grating_cycle,
                swimming_threshold=-2))

        return stimuli


if __name__ == "__main__":
    s = Stytra(protocol=ClosedLoop1DProt())
