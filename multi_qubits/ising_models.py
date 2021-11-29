import numpy as np
import qutip as qt
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import time
import ising_schedules

# Make a base class for doing all common operations on an annealing gadget
class AnnealingGadget:
    def drawHpGraph(self):
        """ Draws the problem Hamiltonian graph.
        """
        fig, ax = plt.subplots(1,1,constrained_layout=True,figsize=(7,7))
        nx.draw(
            self.Hp_graph,
            with_labels=True,
            pos=self.Hp_graph_layout,
            node_color=self.Hp_node_colors,
            node_size=3000,
            width=2,
            labels=self.Hp_node_labels,
            font_weight="bold",
            font_size=16,
            ax=ax
        )
        edge_labels = nx.get_edge_attributes(self.Hp_graph,'desc')
        nx.draw_networkx_edge_labels(
            self.Hp_graph,
            self.Hp_graph_layout,
            edge_labels=edge_labels,
            font_weight="bold",
            font_size=16,
            ax=ax
        )
    
    def drawHdGraph(self):
        """ Draws the driver Hamiltonian graph.
        """
        fig, ax = plt.subplots(1,1,constrained_layout=True,figsize=(7,7))
        nx.draw(
            self.Hd_graph,
            with_labels=True,
            pos=self.Hd_graph_layout,
            node_color=self.Hd_node_colors,
            node_size=3000,
            width=2,
            labels=self.Hd_node_labels,
            font_weight="bold",
            font_size=16,
            ax=ax
        )
        edge_labels = nx.get_edge_attributes(self.Hd_graph,'desc')
        nx.draw_networkx_edge_labels(
            self.Hd_graph,
            self.Hd_graph_layout,
            edge_labels=edge_labels,
            font_weight="bold",
            font_size=16,
            ax=ax
        )
    
    def drawHGraphs(self):
        fig, ax = plt.subplots(1,2,constrained_layout=True,figsize=(15,7))
        #fig.suptitle("Driver and Problem Graphs")
        
        nx.draw(
            self.Hd_graph,
            with_labels=True,
            pos=self.Hd_graph_layout,
            node_color=self.Hd_node_colors,
            node_size=3000,
            width=2,
            labels=self.Hd_node_labels,
            font_weight="bold",
            font_size=16,
            ax=ax[0]
        )
        edge_labels = nx.get_edge_attributes(self.Hd_graph,'desc')
        nx.draw_networkx_edge_labels(
            self.Hd_graph,
            self.Hd_graph_layout,
            edge_labels=edge_labels,
            font_weight="bold",
            font_size=16,
            ax=ax[0]
        )
        ax[0].set_title("Driver Graph")
        
        nx.draw(
            self.Hp_graph,
            with_labels=True,
            pos=self.Hp_graph_layout,
            node_color=self.Hp_node_colors,
            node_size=3000,
            width=2,
            labels=self.Hp_node_labels,
            font_weight="bold",
            font_size=16,
            ax=ax[1]
        )
        edge_labels = nx.get_edge_attributes(self.Hp_graph,'desc')
        nx.draw_networkx_edge_labels(
            self.Hp_graph,
            self.Hp_graph_layout,
            edge_labels=edge_labels,
            font_weight="bold",
            font_size=16,
            ax=ax[1]
        )
        ax[1].set_title("Problem Graph")
    
    def annealInstantaneousGap(self, sparams, tan_pts, gaps=None, timeit=False, rev=False):
        """ Runs an annealing schedule and gets the instantaneous transition energies.
        """
        if timeit:
            start = time.time()
        if gaps is None:
            gaps = self.Hd.shape[0]-1 # Get all gaps
        
        Egaps = np.array([[0.0]*len(tan_pts)]*(gaps+1))
        for i,t in enumerate(tan_pts):
            Ei = (ising_schedules.A(t, sparams)*self.Hd + ising_schedules.B(t, sparams)*self.Hp).eigenenergies() if not rev else (ising_schedules.A(t, sparams)*self.Hd + ising_schedules.B(t, sparams)*self.Hp).eigenenergies()[::-1]
            for j in range(gaps+1):
                Egaps[j,i] = Ei[j]-Ei[0]
        if timeit:
            duration = time.time()-start
            print ("Duration:\t%.1f s" % duration)
        return Egaps

    def annealMinimumGaps(self, tan_pts, Egaps):
        """ Finds where the minimum instantaneous transition energies occur in normalised time.
        """
        minima = [0]*(len(Egaps)-1)
        for i,Epts in enumerate(Egaps):
            # Discard first point (no gap)
            if i == 0:
                continue
            
            # Find minimum index and save xy point
            j = np.argmin(Epts)
            minima[i-1] = [tan_pts[j],Epts[j]]
        return np.array(minima)
    
    def annealInstantaneousOverlaps(self, sparams, tan_pts, timeit=False, coefs=5, signed=False):
        """ Gets coefficients of each eigenstate in terms of the problem solution state
        """
        if timeit:
            start = time.time()
        
        # Get problem solution
        Ep,Vp = self.Hp.eigenstates()
        
        # Get anneal states
        ci2 = np.array([[0.0]*len(tan_pts)]*coefs)
        Ei = None
        Vi = None
        for i,t in enumerate(tan_pts):
            Ei,Vi = (ising_schedules.A(t, sparams)*self.Hd + ising_schedules.B(t, sparams)*self.Hp).eigenstates()
            
            # Project
            for j,Vik in enumerate(Vi[:coefs]):
                ci = (Vp[0].dag()*Vik).data.todense()[0,0]
                ci2[j,i] = np.sign(ci)*np.absolute(ci)**2 if signed else np.absolute(ci)**2
        if timeit:
            duration = time.time()-start
            print ("Duration:\t%.1f s" % duration)
        return ci2
    
    def annealInstantaneousPopulations(self, sparams, tan_pts, timeit=False, coefs=5, signed=False):
        """ Gets coefficients of each eigenstate in terms of the problem solution state
        """
        if timeit:
            start = time.time()
        
        # Get problem solution
        Ep,Vp = self.Hp.eigenstates()
        
        # Get anneal states
        ci2 = np.array([[0.0]*len(tan_pts)]*coefs)
        Ei = None
        Vi = None
        for i,t in enumerate(tan_pts):
            Ei,Vi = (ising_schedules.A(t, sparams)*self.Hd + ising_schedules.B(t, sparams)*self.Hp).eigenstates()
            
            # Project
            for j,Vpk in enumerate(Vp[:coefs]):
                ci = (Vi[0].dag()*Vpk).data.todense()[0,0]
                ci2[j,i] = np.sign(ci)*np.absolute(ci)**2 if signed else np.absolute(ci)**2
        if timeit:
            duration = time.time()-start
            print ("Duration:\t%.1f s" % duration)
        return ci2
    
    def annealStateProbability(self, sparams, tan_pts, states=None, timeit=False):
        """ Finds the probabilities of finding the states of the system as a function of total annealing time.
        """
        if timeit:
            start = time.time()
        if states is None:
            states = self.Hd.shape[0] # Get all states
        
        # Initial Hamiltonian
        H0 = ising_schedules.A(0, sparams)*self.Hd + ising_schedules.B(0, sparams)*self.Hp
        
        # Initialise in the ground state
        E,V = H0.eigenstates()
        Einit = E[0]
        Vinit = V[0]
        
        # Use problem Hamiltonian ground state
        Ep,Vp = self.Hp.eigenstates()
        
        # Annealing H
        H = [[self.Hd,ising_schedules.A],[self.Hp,ising_schedules.B]]
        
        # Solve the time dependence
        p = np.array([np.zeros(len(tan_pts))]*states)
        for i,tan in enumerate(tan_pts):
            times = np.array([0,tan]) # Solve for final time only, initial time must be specified too
            
            # Use time dependent Hamiltonian and specify initial state
            # Unitary evolution: No dissipation. qt.sigmaz()
            sparams["tan"] = tan
            result = qt.mesolve(H, Vinit, times, [], [], args=sparams, options=qt.Options(nsteps=1000000000))
            
            # Save relevant results
            for j in range(states):
                p[j,i] = qt.expect(qt.ket2dm(result.states[-1]),Vp[j])
        if timeit:
            duration = time.time()-start
            print ("Duration:\t%.1f s" % duration)
        return p

    def annealStateEvolution(self, sparams, t_pts, states=None, timeit=False):
        """ Finds the evolution of the states during an anneal.
        """
        if timeit:
            start = time.time()
        if states is None:
            states = self.Hd.shape[0] # Get all states
        
        # Initial Hamiltonian
        H0 = ising_schedules.A(t_pts[0], sparams)*self.Hd + ising_schedules.B(t_pts[0], sparams)*self.Hp
        
        # Initialise in the ground state
        E,V = H0.eigenstates()
        Einit = E[0]
        Vinit = V[0]
        
        # Use problem Hamiltonian ground state
        Ep,Vp = self.Hp.eigenstates()
        
        # Annealing H
        H = [[self.Hd,ising_schedules.A],[self.Hp,ising_schedules.B]]
        
        # Solve the time dependence
        p = np.array([np.zeros(len(t_pts))]*states)
        result = qt.mesolve(H, Vinit, t_pts, [], [], args=sparams, options=qt.Options(nsteps=1000000000))
        
        # Save relevant results
        for i in range(len(t_pts)):
            for j in range(states):
                p[j,i] = qt.expect(qt.ket2dm(result.states[i]),Vp[j])
        if timeit:
            duration = time.time()-start
            print ("Duration:\t%.1f s" % duration)
        return p
    
    def initHamiltonian(self):
        """ Initialises the driver and problem Hamiltonians given the associated graphs.
        """
        
        # HS expander
        I2 = qt.qeye(2)
        
        # Get bias parameters
        qubit_hz = nx.get_node_attributes(self.Hp_graph,'hz')
        qubit_hx = nx.get_node_attributes(self.Hd_graph,'hx')
        
        # Setup qubits
        self.Hp = 0
        self.Hd = 0
        for i in range(self.Nqb):
            fullx = [I2]*self.Nqb
            fullz = [I2]*self.Nqb
            
            fullx[i] = qubit_hx[i]*qt.sigmax()
            fullz[i] = qubit_hz[i]*qt.sigmaz()
            
            self.Hd += qt.tensor(*fullx)
            self.Hp += qt.tensor(*fullz)
        
        # Setup couplers
        Hp_edge_atts = nx.get_edge_attributes(self.Hp_graph,'J')
        Hp_edges = list(self.Hp_graph.edges)
        Hd_edge_atts = nx.get_edge_attributes(self.Hd_graph,'J')
        Hd_edges = list(self.Hd_graph.edges)
        
        
        for i in range(len(Hp_edges)):
            fullz = [I2]*self.Nqb
            
            edge = Hp_edges[i]
            fullz[edge[0]] = qt.sigmaz()
            fullz[edge[1]] = qt.sigmaz()
            
            self.Hp += Hp_edge_atts[edge]*qt.tensor(*fullz)
        
        for i in range(len(Hd_edges)):
            fullx = [I2]*self.Nqb
            
            edge = Hd_edges[i]
            fullx[edge[0]] = qt.sigmay()
            fullx[edge[1]] = qt.sigmay()
            
            self.Hd += Hd_edge_atts[edge]*qt.tensor(*fullx)
        
        # Convert to angular frequency
        self.Hp *= 2*np.pi
        self.Hd *= 2*np.pi
    
    def getHp(self):
        return self.Hp

    def getHd(self):
        return self.Hd

# Tameem Triangle Gadget
class TriangleGadget(AnnealingGadget):
    pass

# Tameem Chain Gadget
class ChainGadget(AnnealingGadget):
    def __init__(self, Nqb, R=1, D=0.05, hx_init=5):
        self.Nqb = Nqb
        
        # Problem graph
        self.Hp_graph = nx.Graph()
        self.Hp_node_colors = []
        self.Hp_node_labels = {}
        
        # Driver graph
        self.Hd_graph = nx.Graph()
        self.Hd_node_colors = []
        self.Hd_node_labels = {}
        
        # Create problem graph
        for qb in range(Nqb//2):
            qb1 = qb*2
            qb2 = qb*2+1
            self.Hp_graph.add_node(qb1,hz=(R*(1-D)))
            self.Hp_graph.add_node(qb2,hz=-R)
            self.Hp_node_colors.append('C0')
            self.Hp_node_colors.append('C1')
            self.Hp_node_labels[qb1] = '$R(1-D)$'
            self.Hp_node_labels[qb2] = '$-R$'
            self.Hp_graph.add_edge(qb1,qb2,J=-R,desc="$-R$")
        for e in range(0,Nqb-1,2):
            if e < Nqb-3:
                self.Hp_graph.add_edge(e,e+2,J=-R,desc="$-R$")
        
        # Create driver graph
        for qb in range(Nqb//2):
            qb1 = qb*2
            qb2 = qb*2+1
            self.Hd_graph.add_node(qb1,hx=hx_init)
            self.Hd_graph.add_node(qb2,hx=hx_init)
            self.Hd_node_colors.append('C0')
            self.Hd_node_colors.append('C1')
            self.Hd_node_labels[qb1] = '$h_x$'
            self.Hd_node_labels[qb2] = '$h_x$'
        
        # Set the graph layout
        self.Hp_graph_layout = nx.circular_layout(self.Hp_graph)
        self.Hd_graph_layout = nx.circular_layout(self.Hd_graph)
        
        # Initialise the driver and problem Hamiltonians
        self.initHamiltonian()
    
# Tameem Loop Gadget
class LoopGadget(AnnealingGadget):
    def __init__(self, Nqb, R=1, hx_init=5, biased_qbs=2):
        self.Nqb = Nqb
        
        # Problem graph
        self.Hp_graph = nx.Graph()
        self.Hp_node_colors = []
        self.Hp_node_labels = {}
        
        # Driver graph
        self.Hd_graph = nx.Graph()
        self.Hd_node_colors = []
        self.Hd_node_labels = {}
        
        # Create problem graph
        for qb in range(Nqb):
            self.Hp_graph.add_node(qb,hz=0)
            self.Hp_node_colors.append('C0')
            self.Hp_node_labels[qb] = '$0$'
        for e in range(Nqb):
            if e == Nqb-1:
                self.Hp_graph.add_edge(0,e,J=-R,desc="$-R$")
            else:
                self.Hp_graph.add_edge(e,e+1,J=-R,desc="$-R$")
        
        # Create driver graph
        for qb in range(Nqb):
            self.Hd_graph.add_node(qb,hx=hx_init)
            self.Hd_node_colors.append('C0')
            self.Hd_node_labels[qb] = '$h_x$'
        
        # Select required number of biased qubits
        # 2 is the only current case
        if biased_qbs == 2:
            self.Hp_graph.nodes[0]['hz'] = -R
            self.Hp_node_colors[0] = 'C1'
            self.Hd_node_colors[0] = 'C1'
            self.Hp_node_labels[0] = '$-R$'
            self.Hp_graph.nodes[self.Nqb//2]['hz'] = R - 1
            self.Hp_node_colors[self.Nqb//2] = 'C2'
            self.Hd_node_colors[self.Nqb//2] = 'C2'
            self.Hp_node_labels[self.Nqb//2] = '$R-1$'
        else:
            raise Exception("can only bias two qubits currently.")
        
        # Set the graph layout
        self.Hp_graph_layout = nx.circular_layout(self.Hp_graph)
        self.Hd_graph_layout = nx.circular_layout(self.Hd_graph)
        
        # Init the coupler strengths near selected biased qubits
        left_edges = self.Hp_graph.edges([0])
        d1 = {}
        d2 = {}
        for e in left_edges:
            d1[e] = -R/2
            d2[e] = "$-R/2$"
        nx.set_edge_attributes(self.Hp_graph,d1,'J')
        nx.set_edge_attributes(self.Hp_graph,d2,'desc')
        
        # Initialise the driver and problem Hamiltonians
        self.initHamiltonian()

# Weighted MIS Problem with XX driver option
# Vicky Choi
class WeightedMIS(AnnealingGadget):
    
    def __init__(self, Hd_node_dict, Hd_edge_dict, Hp_node_dict, Hp_edge_dict):
        
        self.Nqb = len(Hp_node_dict)
        
        # Problem graph
        self.Hp_graph = nx.Graph()
        self.Hp_node_colors = []
        self.Hp_node_labels = {}
        
        # Driver graph
        self.Hd_graph = nx.Graph()
        self.Hd_node_colors = []
        self.Hd_node_labels = {}
        
        # Create problem graph
        for n, hz in Hp_node_dict.items():
            self.Hp_graph.add_node(n,hz=hz)
            self.Hp_node_colors.append('C0')
            self.Hp_node_labels[n] = '$%.3f$'%(hz)
        for e, Jzz in Hp_edge_dict.items():
            self.Hp_graph.add_edge(e[0],e[1],J=Jzz,desc=('$%.3f$'%Jzz))
        
        # Create driver graph
        for n, hx in Hd_node_dict.items():
            self.Hd_graph.add_node(n,hx=hx)
            self.Hd_node_colors.append('C0')
            self.Hd_node_labels[n] = '$%.3f$'%(hx)
        for e, Jxx in Hd_edge_dict.items():
            self.Hd_graph.add_edge(e[0],e[1],J=Jxx,desc=('$%.3f$'%Jxx))
        
        # Set the graph layout
        self.Hp_graph_layout = nx.spring_layout(self.Hp_graph)
        self.Hd_graph_layout = nx.spring_layout(self.Hd_graph)
        
        # Initialise the driver and problem Hamiltonians
        self.initHamiltonian()
        #self.Hp =  self.Hp
        #self.Hd =  self.Hd
    
    
