## Single Qubit Annealing Notebooks

A collection of Jupyter notebooks that explore Single Qubit Annealing using python library `qutip` for defining Hamiltonians and solving time dependent solutions.

#### Available Files:

 * `sqb_oscillations.ipynb`
   - Annealing of spin-level Hamiltonians.
   - Exploration of different annealing schedules.
   - Use of the WKB approximation for approximating non-linearity of tunneling term.
   - Exploration of bias line circuit effects.
   - Annealing a CSFQ circuit Hamiltonian.
   
**TODO:** Keep this notebook specific to spin-level Hamiltonians. Circuit effects should be included in circuit model of qubits.

**TODO:** Add decoherence models (probably not necessary at this stage, maybe save it for more advanced models).
 
 * `LZS_interference.ipynb`
   - Spin-level investigation of Landau-Zener-Stückelberg interference effects.
   - Comparison with half-wave schedules.

**TODO:** Synchronise with later version where errors and mistakes were corrected.

**TODO:** Determine the probabilities of LZ tunneling events.

 * `csfq_annealing.ipynb`
   - CSFQ circuit model investigation of single qubit annealing.

**TODO:** Decentralise class and function definitions from this notebook and modularise.

**TODO:** Include bias-line circuit effects.

**TODO:** Use accurate mapping technique (SW is appealing).

**TODO:** Implement Hamiltonian time evolution solver for full circuit model (small time steps).

 * `csfq_resonator_interaction.ipynb`
   - CSFQ-resonator capacitively coupled model investigation.
   - Use of Mostafa's JC ladder code to solve resonator shifts.
   - Development of data fitting procedure for determining system parameters from resonator shifts.

**TODO:** Decentralise class and function definitions from this notebook and modularise.

**TODO:** Finish data fitting procedure.

Ideas:

 - Lupascu: Adaptive annealing when approaching minimum gap using weak measurement
 - The anneal is slowed down when approaching the minimum gap to suppress diabatic transitions.
 - Can this be used in a different way to speed up the anneal?
   * Do a first anneal to detect where the gap is and what it looks like.
   * Then do a second anneal where the schedules are modified to turn the gap into a double gap.
   * Then with knowledge of the shape and location of the gap, deliberatly create a double LZ event and exploit interference.
