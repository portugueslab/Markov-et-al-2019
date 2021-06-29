import numpy as np
import pandas as pd

from stytra_config import ConfiguredStytra
from stytra.stimulation import Protocol
from stytra.stimulation.stimuli import CalibratingClosedLoop1D, \
    GratingStimulus, Basic_CL_1D, GainLagClosedLoop1D, GainChangerStimulus
from lightparam import Param


class ClosedLoop1DProt(Protocol):
    name = "E0030_long_term_adaptation/v09_lightsheet"

    stytra_config = dict(
        tracking=dict(embedded=True, method="tail", estimator="vigor"),
        log_format="hdf5"
        )

    def __init__(self):
        super().__init__()

        self.grating_cycle = Param(10, loadable=False, editable=False)
        self.n_repeats_selfcalib = Param(10, loadable=False, editable=False)
        self.n_repeats_pre = Param(10, loadable=False, editable=False)
        self.n_repeats_exp = Param(50, loadable=False, editable=False)
        self.n_repeats_post = Param(50, loadable=False, editable=False)
        self.target_vel = Param(-20., limits=(-20, 20), loadable=False, editable=False)
        self.calibrate_after = Param(10, loadable=False, editable=False)
        self.max_interbout_time = Param(300, loadable=False, editable=False)
        self.min_bouts_n = Param(20, loadable=False, editable=False)
        self.start_gain = Param(-26, limits=(-50., 0.), loadable=False, editable=False)
        self.lag = Param(0.225, limits=[0, 0.225])
        self.swimming_threshold = Param(-5, limits=(-20, 20))
        inter_stim_pause = 7.5
        grating_vel = 10
        grating_duration = 15
        self.inter_stim_pause = Param(inter_stim_pause, loadable=False, editable=False)
        self.grating_vel = Param(grating_vel, loadable=False, editable=False)
        self.grating_duration = Param(grating_duration, loadable=False, editable=False)

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

        ClosedLoop1DGratings = type("Stim", (GainLagClosedLoop1D,
                                             GratingStimulus), {})
        lag = float(self.lag)
        print(lag)
        return ClosedLoop1DGratings(
                df_param=df,
                lag=lag,
                grating_angle=np.pi / 2,
                grating_period=self.grating_cycle,
                swimming_threshold=self.swimming_threshold,
                max_interbout_time=self.max_interbout_time)

    def get_stim_sequence(self):
        stimuli = []
        # # gratings
        # Calibration:
        stimuli.append(GainChangerStimulus(self.start_gain))

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
                min_bout_n=self.min_bouts_n,
                swimming_threshold=self.swimming_threshold,
                target_avg_fish_vel=self.target_vel,
                calibrate_after=self.calibrate_after,
                max_interbout_time=self.max_interbout_time))

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
                swimming_threshold=self.swimming_threshold,
                max_interbout_time=self.max_interbout_time))

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
                swimming_threshold=self.swimming_threshold,
                max_interbout_time=self.max_interbout_time))

        return stimuli


if __name__ == "__main__":
    s = ConfiguredStytra(protocol=ClosedLoop1DProt())
