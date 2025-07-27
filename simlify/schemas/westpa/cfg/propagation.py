# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scientific Computing Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

from pydantic import BaseModel, Field


class PropagationConfig(BaseModel):
    """Propagation settings for WESTPA."""

    gen_istates: bool = Field(default=True)
    """Boolean specifying whether to generate initial states from the basis states.
    The executable propagator defines a specific configuration block, and custom
    propagators should override the `WESTPropagator.gen_istate()` method.
    """

    block_size: int = Field(default=1)
    """An integer defining how many segments should be passed to a worker at a time.
    When using the serial work manager, this value should be set to the maximum number
    of segments per iteration to avoid significant overhead incurred by the locking
    mechanism in the WMFutures framework. Parallel work managers might benefit from
    setting this value greater than one in some instances to decrease network
    communication load.
    """
    save_transition_matrices: bool = True
    """Save transition matrices."""
    max_run_wallclock: str | None = None
    """A time in dd:hh:mm:ss or hh:mm:ss specifying the maximum wallclock time of a
    particular WESTPA run. If running on a batch queuing system, this time should be
    set to less than the job allocation time to ensure that WESTPA shuts down cleanly.
    """
    max_total_iterations: int | None = None
    """Max number of iterations to run."""
