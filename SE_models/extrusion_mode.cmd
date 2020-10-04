// extrusion_mode.fe
//
// This is a command file for extru_geom.fe
// Use Surface Evolver (http://facstaff.susqu.edu/brakke/evolver/evolver.html) to call extru_geom.fe
//
// Surface Evolver Commands
// - generate geometry of the center cell
// - generate geometry of the ring of neighbours
//
// By: Lester Lin (lesterlin88@gmail.com)
//****************************************************************************************************

// set the discretization pattern
// 'api_bnumber' specifies the number of apical polar bands to build
api_bnumber:= 3;
z_height := 20;

// set the discretization pattern
// this is to specify the number of sides for the cell
discret_max := 18;
discret_half := discret_max/2;

//****************************************************************************************************
//********************
// Build center cell
//********************

// Build basal
buildbasal := {

	/*
	// to build the basal center vertex
	new_basalvertex := new_vertex(0,
	                           	   0,
							       0);
	set vertex[new_basalvertex] base_center 1;
	set vertex[new_basalvertex] constraint 2;		*/

	// Use a non-negative constraint here so that the basal vertex point can be lifted up during
	// 2nd phase of extrusion. If base pole vertex is constrained at z=0 at first, then I need to
	// unset the constraint 1 before going on to 2nd phase.

	// Form the basal ring first, which is also the basal ring to be colored in magenta
	// create new outermost basal vertices
	for ( bb := 1 ; bb < (discret_max+1) ; bb := bb + 1 ) {
		new_basevertex := new_vertex((major)*cos(bb*pi/discret_half),
									 (minor)*sin(bb*pi/discret_half),
									 0
									 );
		set vertex[new_basevertex] v_side bb;
		set vertex[new_basevertex] cringvertex 1;
		set vertex[new_basevertex] constraint 1;
	};

	// Loop to build edges only for the initial cringedge==1
	for ( bb := 1 ; bb < (discret_max+1) ; bb := bb + 1 ) {

		if (bb == discret_max) then {
			foreach vertex vv where (cringvertex == 1 && v_side == bb) do {
				foreach vertex vvv where (cringvertex == 1 && v_side == 1) do {
					new_basering := new_edge(vv.id, vvv.id);
					set edge[new_basering] cringedge 1;
					set edge[new_basering] e_side bb;
				};
			};
		}
		else {
			foreach vertex vv where (cringvertex == 1 && v_side == bb) do {
				foreach vertex vvv where (cringvertex == 1 && v_side == (bb+1)) do {
					new_basering := new_edge(vv.id, vvv.id);
					set edge[new_basering] cringedge 1;
					set edge[new_basering] e_side bb;
				};
			};
		};
	};

	// this step is to first pick out the basal ring using the 'magenta' color, so that it can be found easily afterwards
	foreach edge eon where (eon.cringedge == 1) do {
		set edge[eon.id] color magenta;
		set edge[eon.id] tension cbasalringtens;		// the tension of basal ring is an initial condition to decide on
	};

for ( bubase := 2 ; bubase < major ; bubase := bubase + 2 ) {

	dissolve facets where (basespoke_tri == 1);
	dissolve edges where (base_spoke == 1);

	// to increment the indexes of attributes above & including the original basal ring (cringvertex == 1)
	foreach vertex vv where (vv.cringvertex > 0) do {
		set vertex[vv.id] cringvertex (vv.cringvertex + 1);
	};
	foreach edge eggy where (eggy.cbandwidth > 0) do {
		set edge[eggy.id] cbandwidth (eggy.cbandwidth + 1);
	};
	foreach edge eggy where (eggy.cringedge > 0) do {
		set edge[eggy.id] cringedge (eggy.cringedge + 1);
	};
	foreach edge eggy where (eggy.bandsq > 0) do {
		set edge[eggy.id] bandsq (eggy.bandsq + 1);
	};
	foreach facet ff where (ff.tri > 0) do {
		set facet[ff.id] tri (ff.tri + 1);
	};

	// create new basal vertices
	for ( bb := 1 ; bb < (discret_max+1) ; bb := bb + 1 ) {
		new_basevertex := new_vertex((major-bubase)*cos(bb*pi/discret_half),
									 (minor-bubase*(minor/major))*sin(bb*pi/discret_half),
									 0
									 );
		set vertex[new_basevertex] v_side bb;
		set vertex[new_basevertex] cringvertex 1;
		set vertex[new_basevertex] constraint 1;
	};

	// Loop to build edges
	for ( bb := 1 ; bb < (discret_max+1) ; bb := bb + 1 ) {

	if (bb == discret_max) then {
		foreach vertex vv where (cringvertex == 1 && v_side == bb) do {
			foreach vertex vvv where (cringvertex == 2 && v_side == bb) do {
				new_basevertical := new_edge(vv.id, vvv.id);
				set edge[new_basevertical] cbandwidth 1;
				set edge[new_basevertical] e_side bb;
			};
		};

		foreach vertex vv where (base_center == 1) do {
			foreach vertex vvv where (cringvertex == 1 && v_side == bb) do {
				new_spoke := new_edge(vv.id, vvv.id);
				set edge[new_spoke] cbandwidth 0;
				set edge[new_spoke] base_spoke 1;
				set edge[new_spoke] e_side bb;
				};
		};

		foreach vertex vv where (cringvertex == 1 && v_side == bb) do {
			foreach vertex vvv where (cringvertex == 1 && v_side == 1) do {
				new_basering := new_edge(vv.id, vvv.id);
				set edge[new_basering] cringedge 1;
				set edge[new_basering] e_side bb;
				};
		};

		foreach vertex vv where (cringvertex == 2 && v_side == 1) do {
			foreach vertex vvv where (cringvertex == 1 && v_side == bb) do {
				new_basediag := new_edge(vv.id, vvv.id);
				set edge[new_basediag] e_diag 1;
				set edge[new_basediag] e_side bb;
				set edge[new_basediag] bandsq 1;
				};
		};
	}
	else {
		foreach vertex vv where (cringvertex == 1 && v_side == bb) do {
			foreach vertex vvv where (cringvertex == 2 && v_side == bb) do {
				new_basevertical := new_edge(vv.id, vvv.id);
				set edge[new_basevertical] cbandwidth 1;
				set edge[new_basevertical] e_side bb;
			};
		};

		foreach vertex vv where (base_center == 1) do {
			foreach vertex vvv where (cringvertex == 1 && v_side == bb) do {
				new_spoke := new_edge(vv.id, vvv.id);
				set edge[new_spoke] cbandwidth 0;
				set edge[new_spoke] base_spoke 1;
				set edge[new_spoke] e_side bb;
				};
		};

		foreach vertex vv where (cringvertex == 1 && v_side == bb) do {
			foreach vertex vvv where (cringvertex == 1 && v_side == (bb+1)) do {
				new_basering := new_edge(vv.id, vvv.id);
				set edge[new_basering] cringedge 1;
				set edge[new_basering] e_side bb;
				};
		};

		foreach vertex vv where (cringvertex == 2 && v_side == (bb+1)) do {
			foreach vertex vvv where (cringvertex == 1 && v_side == bb) do {
				new_basediag := new_edge(vv.id, vvv.id);
				set edge[new_basediag] e_diag 1;
				set edge[new_basediag] e_side bb;
				set edge[new_basediag] bandsq 1;
				};
		};
	};
	};

	// Loop to build facets
	for ( bb := 1 ; bb < (discret_max+1) ; bb := bb + 1 ) {

	if (bb == discret_max) then {
		// Right triangle
		foreach edge ee1 where (cbandwidth == 1 && e_side == 1) do {
			foreach edge ee2 where (bandsq == 1 && e_side == bb) do {
				foreach edge ee3 where (cringedge == 1 && e_side == bb) do {
					new_right_tri := new_facet(ee1.id, ee2.id, ee3.id);
					set facet[new_right_tri] tri 1;
					set facet[new_right_tri] f_side bb;
				};
			};
		};

		// Left triangle
		foreach edge ee1 where (bandsq == 1 && e_side == bb) do {
			foreach edge ee2 where (cringedge == 2 && e_side == bb) do {
				foreach edge ee3 where (cbandwidth == 1 && e_side == bb) do {
					new_left_tri := new_facet(-ee1.id, -ee2.id, -ee3.id);
					set facet[new_left_tri] tri 1;
					set facet[new_left_tri] f_side bb;
				};
			};
		};

		// Spoke triangle
		foreach edge ee1 where (base_spoke == 1 && e_side == 1) do {
			foreach edge ee2 where (cringedge == 1 && e_side == bb) do {
				foreach edge ee3 where (base_spoke == 1 && e_side == bb) do {
					new_spoke_tri := new_facet(ee1.id, -ee2.id, -ee3.id);
					set facet[new_spoke_tri] tri 0;
					set facet[new_spoke_tri] basespoke_tri 1;
					set facet[new_spoke_tri] f_side bb;
				};
			};
		};
	}
	else {
		// Right triangle
		foreach edge ee1 where (cbandwidth == 1 && e_side == (bb+1)) do {
			foreach edge ee2 where (bandsq == 1 && e_side == bb) do {
				foreach edge ee3 where (cringedge == 1 && e_side == bb) do {
					new_right_tri := new_facet(ee1.id, ee2.id, ee3.id);
					set facet[new_right_tri] tri 1;
					set facet[new_right_tri] f_side bb;
				};
			};
		};

		// Left triangle
		foreach edge ee1 where (bandsq == 1 && e_side == bb) do {
			foreach edge ee2 where (cringedge == 2 && e_side == bb) do {
				foreach edge ee3 where (cbandwidth == 1 && e_side == bb) do {
					new_left_tri := new_facet(-ee1.id, -ee2.id, -ee3.id);
					set facet[new_left_tri] tri 1;
					set facet[new_left_tri] f_side bb;
				};
			};
		};

		// Spoke triangle
		foreach edge ee1 where (base_spoke == 1 && e_side == (bb+1)) do {
			foreach edge ee2 where (cringedge == 1 && e_side == bb) do {
				foreach edge ee3 where (base_spoke == 1 && e_side == bb) do {
					new_spoke_tri := new_facet(ee1.id, -ee2.id, -ee3.id);
					set facet[new_spoke_tri] tri 0;
					set facet[new_spoke_tri] basespoke_tri 1;
					set facet[new_spoke_tri] f_side bb;
				};
			};
		};
	}; // end of else statement for building facets
	}; // end of loop for building facets

}; // end of 'bubase' for-loop

}; // end of 'buildbasal' command
buildbasal;

//****************************************************************************************************
// Build lateral
buildlat := {

	// Note: use 'bulat' to determine the height and z-interval of each band
	for ( bulat := 2 ; bulat <= 20 ; bulat := bulat + 2 ){

		// to count the current total number of bands excluding the basal pole
		band_max := max(edge elsa where (elsa.cbandwidth > 0), elsa.cbandwidth);
		ring_max := max(edge elsa where (elsa.cringedge > 0), elsa.cringedge);
		tri_max := max(facet ff where (ff.tri > 0), ff.tri);

		// create new lateral vertices
		for ( bb := 1 ; bb < (discret_max+1) ; bb := bb + 1 ) {
		new_cringvertex := new_vertex((major)*cos(bb*pi/discret_half),
									 (minor)*sin(bb*pi/discret_half),
									 (bulat)
									 );
		set vertex[new_cringvertex] v_side bb;
		set vertex[new_cringvertex] cringvertex (ring_max + 1);
	};

		// Loop to build edges
		for ( bb := 1 ; bb < (discret_max+1) ; bb := bb + 1 ) {

		if (bb == discret_max) then {
			// vertical
			foreach vertex vv where (cringvertex == (ring_max) && v_side == bb) do {
				foreach vertex vvv where (cringvertex == (ring_max + 1) && v_side == bb) do {
					new_vertical := new_edge(vv.id, vvv.id);
					set edge[new_vertical] cbandwidth (band_max + 1);
					set edge[new_vertical] e_side bb;
				};
			};
			// diagonal edge for new band, pointing from top right to bottom left
			foreach vertex vv where (cringvertex == (ring_max + 1) && v_side == 1) do {
				foreach vertex vvv where (cringvertex == ring_max && v_side == bb) do {
					new_diag := new_edge(vv.id, vvv.id);
					set edge[new_diag] e_diag 1;
					set edge[new_diag] e_side bb;
					set edge[new_diag] bandsq (band_max + 1);
				};
			};
			// new cringedge
			foreach vertex vv where (cringvertex == (ring_max + 1) && v_side == bb) do {
				foreach vertex vvv where (cringvertex == (ring_max + 1) && v_side == 1) do {
					new_horizontal := new_edge(vv.id, vvv.id);
					set edge[new_horizontal] cringedge (ring_max + 1);
					set edge[new_horizontal] e_side bb;
				};
			};
		}
		else {
			// vertical
			foreach vertex vv where (cringvertex == (ring_max) && v_side == bb) do {
				foreach vertex vvv where (cringvertex == (ring_max + 1) && v_side == bb) do {
					new_vertical := new_edge(vv.id, vvv.id);
					set edge[new_vertical] cbandwidth (band_max + 1);
					set edge[new_vertical] e_side bb;
				};
			};
			// diagonal edge for new band, pointing from top right to bottom left
			foreach vertex vv where (cringvertex == (ring_max + 1) && v_side == (bb+1)) do {
				foreach vertex vvv where (cringvertex == ring_max && v_side == bb) do {
					new_diag := new_edge(vv.id, vvv.id);
					set edge[new_diag] e_diag 1;
					set edge[new_diag] e_side bb;
					set edge[new_diag] bandsq (band_max + 1);
				};
			};
			// new cringedge
			foreach vertex vv where (cringvertex == (ring_max + 1) && v_side == bb) do {
				foreach vertex vvv where (cringvertex == (ring_max + 1) && v_side == (bb+1)) do {
					new_horizontal := new_edge(vv.id, vvv.id);
					set edge[new_horizontal] cringedge (ring_max + 1);
					set edge[new_horizontal] e_side bb;
				};
			};
		}; // end of else for edge
	}; // end of building edges

		// Loop to build facets
		for ( bb := 1 ; bb < (discret_max+1) ; bb := bb + 1 ) {

		if (bb == discret_max) then {
		// Right triangle on new band
		foreach edge ee1 where (cbandwidth == (band_max + 1) && e_side == 1) do {
			foreach edge ee2 where (bandsq == (band_max + 1) && e_side == bb) do {
				foreach edge ee3 where (cringedge == (ring_max) && e_side == bb) do {
					new_right_tri := new_facet(ee1.id, ee2.id, ee3.id);
					set facet[new_right_tri] tri (band_max + 1);
					set facet[new_right_tri] f_side bb;
				};
			};
		};
		// Left triangle on new band
		foreach edge ee1 where (bandsq == (band_max + 1) && e_side == bb) do {
			foreach edge ee2 where (cringedge == (ring_max + 1) && e_side == bb) do {
				foreach edge ee3 where (cbandwidth == (band_max + 1) && e_side == bb) do {
					new_left_tri := new_facet(-ee1.id, -ee2.id, -ee3.id);
					set facet[new_left_tri] tri (band_max + 1);
					set facet[new_left_tri] f_side bb;
				};
			};
		};
		}
		else {
		// Right triangle on new band
		foreach edge ee1 where (cbandwidth == (band_max + 1) && e_side == (bb+1)) do {
			foreach edge ee2 where (bandsq == (band_max + 1) && e_side == bb) do {
				foreach edge ee3 where (cringedge == (ring_max) && e_side == bb) do {
					new_right_tri := new_facet(ee1.id, ee2.id, ee3.id);
					set facet[new_right_tri] tri (band_max + 1);
					set facet[new_right_tri] f_side bb;
				};
			};
		};
		// Left triangle on new band
		foreach edge ee1 where (bandsq == (band_max + 1) && e_side == bb) do {
			foreach edge ee2 where (cringedge == (ring_max + 1) && e_side == bb) do {
				foreach edge ee3 where (cbandwidth == (band_max + 1) && e_side == bb) do {
					new_left_tri := new_facet(-ee1.id, -ee2.id, -ee3.id);
					set facet[new_left_tri] tri (band_max + 1);
					set facet[new_left_tri] f_side bb;
				};
			};
		};
	}; // end of else statement for building facets
	}; // end of loop for building facets

}; // end of 'bulat' for-loop

	// to count the current total number of bands excluding the basal pole
	band_max := max(edge elsa where (elsa.cbandwidth > 0), elsa.cbandwidth);
	ring_max := max(edge elsa where (elsa.cringedge > 0), elsa.cringedge);
	tri_max := max(facet ff where (ff.tri > 0), ff.tri);

}; // end of 'buildlat' command
buildlat;

//****************************************************************************************************
