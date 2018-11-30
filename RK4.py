from sympy import *

def d_single_f_term(ft, y_dependent=True):

    """ return dts and dys lists for single fterm """

    dt = [  ft[:] + ['t'] ]
    dy = [ ft[:] + ['y'], [] ] 
    return dt, dy


def d_product(term, y_dependent=True):

	# for each f term in terms (connected by multiplication)
	# apply d_single_f_term, replace index, and add list of terms to term list to return

	# product rule
	n = len(term) # number of f multipliers in term
	dt_terms, dy_terms = [], [] #[term[:]]*n, [term[:]]*n
	for f in range(n):

		dt_pop, dy_pop = term[:], term[:]
		del dt_pop[f]
		del dy_pop[f]
		
		dti, dyi = d_single_f_term(term[f], y_dependent=y_dependent)
		dt_terms += [ dti + dt_pop ]
		dy_terms += [ dyi + dy_pop ]

	return dt_terms, dy_terms


def d_terms(terms, y_dependent=True):

	# differentiate all terms

	dt_terms, dy_terms = [], []

	for term in terms:
		dt_terms_i, dy_terms_i = d_product(term, y_dependent=y_dependent)

		dt_terms += dt_terms_i
		dy_terms += dy_terms_i

	return dt_terms, dy_terms

from sympy import *
init_printing()

def get_f_sympy(f_list_rep):

	f = Function("f")
	t,y = symbols("t,y")


	if len(f_list_rep) > 0:


		sorted_partials = sorted(f_list_rep)

		partial_symbols = [symbols(str(partial)) for partial in sorted_partials]

		return Derivative(f(t,y), *partial_symbols)

	else: 
		return f(t,y)



def get_term_sympy(term):

	term_sympy = 1
	for f in term:
		term_sympy =  term_sympy*get_f_sympy(f)
	return term_sympy


def get_f_string(f):

	if len(f) > 0:
		f_str = 'f_{'

		sorted_partials = sorted(f)
		for partial in sorted_partials:
			f_str += partial

		return f_str + '}'

	else: 
		return 'f'

def get_term_string(term):
	term_str = ''

	for f in term:
		term_str += get_f_string(f)

	return term_str

def get_terms_sympy(terms_list):
    
    expression = 0
    for term in terms_list:
        expression += get_term_sympy(term)
        
    return expression


class Node: 

    def __init__(self, terms_list, delta_t, delta_y, multiplier=1): 
        self.terms_list = terms_list
        self.delta_t, self.delta_y = delta_t, delta_y
        self.multiplier = multiplier

        # node children
        self.dt_terms_node = None
        self.dy_terms_node = None
        
        self.sympy_rep = get_terms_sympy(self.terms_list)
        
    def expand_terms(self, delta_t, delta_y,depth, y_dependent=True, depth_offset=0):
        
        dt_terms, dy_terms = d_terms(self.terms_list, y_dependent=y_dependent)
        self.dt_terms_node, self.dy_terms_node = Node(dt_terms,delta_t, delta_y),  Node(dy_terms,delta_t, delta_y)
        
        # apply multipliers
        self.dt_terms_node.apply_multiplier(delta_t*self.multiplier, depth) #depth+1)
        self.dy_terms_node.apply_multiplier(delta_y*self.multiplier, depth) #depth+1)
        
        
    def apply_multiplier(self, delta, depth):
        self.multiplier = delta
        self.sympy_rep = self.sympy_rep*self.multiplier
        
        # apply factorial multiplier
        self.sympy_rep = self.sympy_rep/factorial(depth)

        
def add_order(node, delta_t, delta_y, depth=0, y_dependent=True, depth_offset=0):
    
    # if leaf node found, expand leaf node
    if node.dt_terms_node is None or node.dy_terms_node is None:
        node.expand_terms(delta_t, delta_y, depth+1, y_dependent=y_dependent, depth_offset=depth_offset)
        
    else:
        add_order(node.dt_terms_node, delta_t, delta_y, depth=depth+1, y_dependent=y_dependent, depth_offset=depth_offset)
        add_order(node.dy_terms_node, delta_t, delta_y, depth=depth+1, y_dependent=y_dependent, depth_offset=depth_offset)
        

def get_combined_terms_list(node):
    if node is not None: 

        if ( node.dt_terms_node is None and node.dy_terms_node is None ): 
            return node.terms_list

        else: 
            return node.terms_list + get_combined_terms_list(node.dt_terms_node) + get_combined_terms_list(node.dy_terms_node)

def get_combined_sympy(node):
    if node is not None: 

        if ( node.dt_terms_node is None and node.dy_terms_node is None ): 
            return node.sympy_rep
        else: 
            return node.sympy_rep + get_combined_sympy(node.dt_terms_node) + get_combined_sympy(node.dy_terms_node)
    
def get_taylor_tree_sympy(root):

    return get_combined_sympy(root) 

class TaylorTree:
    
    def __init__(self, order, delta_t=0, delta_y=0, multiplier=1, depth_coeff_offset=0, y_dependent=True):
        
        self.delta_t, self.delta_y = delta_t, delta_y
        self.root = Node([[[]]], self.delta_t, self.delta_y, multiplier=multiplier)
        self.order = order # plus depth coef offset possibly as an option (or manual option)
        
        # depth of tree is == to the order passed
        for ol in range(self.order-1):
            add_order(self.root, self.delta_t, self.delta_y, depth=depth_coeff_offset, y_dependent=y_dependent, depth_offset=depth_coeff_offset)
            
            
    def get_sympy(self):
        h = symbols("h")
        return expand(get_taylor_tree_sympy(self.root)) + Order(h**self.order)
        
# sympy
f = Function("f")
t,y,h = symbols("t,y,h")
yn = symbols("y_n")

# butcher
b1, b2, b3 = 2/9, 3/9, 4/9
c2, a21 = 1/2, 1/2
c3, a31, a32 = 3/4, 0, 3/4

# k1, k2, k3
k1 = h*f(t,y)

tk2 = TaylorTree(3, h*c2, h*a21)
k2 = h*tk2.get_sympy() 

tk3 = TaylorTree(4, h*c3, h*a32)
k3 = expand(h*tk3.get_sympy())

# expand runge kutta method
rk_expansion = expand(yn + b1*k1 + b2*k2 + b3*k3)

# expand taylor series of correct solution
tt = TaylorTree(3, h, h, depth_coeff_offset=1)
y_talyor_expansion = yn + expand(h*tt.get_sympy())

rk_expansion

y_talyor_expansion

rk_expansion  - y_talyor_expansion
