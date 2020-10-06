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
// Build apical lateral, referring to the lateral elements on the same initial apical plane as
// the apical polar cap
buildapi := {
	for ( buapi := 2 ; buapi < major-apimajor ; buapi := buapi + 2 ){

		dissolve facets where (apispoke_tri == 1);
		dissolve edges where (api_spoke == 1);

		band_max := max(edge elsa where (elsa.cbandwidth > 0), elsa.cbandwidth);
		ring_max := max(edge elsa where (elsa.cringedge > 0), elsa.cringedge);
		tri_max := max(facet ff where (ff.tri > 0), ff.tri);

		cell_height := avg(vertex vv where (vv.cringvertex == ring_max), vv.z);
		//	Equation of ellipse, where 'A' is the major axis on the horizontal x-y plane,
		// 'B' is minor axis in the z-plane
		//	Solving for 'y' would find the z-position for each ring of apical vertices
		//	'(20+sqrt((2^2)*(1 - (((major-buapi)^2)/(20^2)))))',
		// the '20' at the front refers to the max cell height at the lateral
		// ((x^2)/(A^2)) + ((y^2)/(B^2)) = 1
		// ((x^2)/(20^2)) + ((y^2)/(2^2)) = 1
		// (y^2) = (2^2)*(1 - ((x^2)/(20^2)))
		//  y = sqrt((2^2)*(1 - ((x^2)/(20^2))))
		//	y = sqrt((2^2)*(1 - (((major-buapi)^2)/(20^2))))

		// build new apical vertices
		for ( bb := 1 ; bb < (discret_max+1) ; bb := bb + 1 ) {
			new_apivertex := new_vertex((major-buapi)*cos(bb*pi/discret_half),
									    (minor-buapi*(minor/major))*sin(bb*pi/discret_half),
									    20
									    //(20+sqrt((3^2)*(1 - (((major-buapi)^2)/(20^2)))))
									    );
			set vertex[new_apivertex] v_side bb;
			set vertex[new_apivertex] cringvertex (ring_max + 1);
		};

		// Loop to build edges
		for ( bb := 1 ; bb < (discret_max+1) ; bb := bb + 1 ) {

		if (bb == discret_max) then {
			foreach vertex vv where (cringvertex == ring_max && v_side == bb) do {
				foreach vertex vvv where (cringvertex == (ring_max + 1) && v_side == bb) do {
					new_apivertical := new_edge(vv.id, vvv.id);
					set edge[new_apivertical] cbandwidth (band_max + 1);
					set edge[new_apivertical] e_side bb;
				};
			};

			foreach vertex vv where (cringvertex == (ring_max + 1) && v_side == bb) do {
				foreach vertex vvv where (api_center == 1) do {
					new_spoke := new_edge(vv.id, vvv.id);
					set edge[new_spoke] cbandwidth (band_max + 2);
					set edge[new_spoke] api_spoke 1;
					set edge[new_spoke] e_side bb;
				};
			};

			foreach vertex vv where (cringvertex == (ring_max + 1) && v_side == bb) do {
				foreach vertex vvv where (cringvertex == (ring_max + 1) && v_side == 1) do {
					new_apiring := new_edge(vv.id, vvv.id);
					set edge[new_apiring] cringedge (ring_max + 1);
					set edge[new_apiring] e_side bb;
				};
			};

			foreach vertex vv where (cringvertex == (ring_max + 1) && v_side == 1) do {
				foreach vertex vvv where (cringvertex == ring_max && v_side == bb) do {
					new_apidiag := new_edge(vv.id, vvv.id);
					set edge[new_apidiag] e_diag 1;
					set edge[new_apidiag] e_side bb;
					set edge[new_apidiag] bandsq (band_max + 1);
				};
			};
		}
		else {
			foreach vertex vv where (cringvertex == ring_max && v_side == bb) do {
				foreach vertex vvv where (cringvertex == (ring_max + 1) && v_side == bb) do {
					new_apivertical := new_edge(vv.id, vvv.id);
					set edge[new_apivertical] cbandwidth (band_max + 1);
					set edge[new_apivertical] e_side bb;
				};
			};

			foreach vertex vv where (cringvertex == (ring_max + 1) && v_side == bb) do {
				foreach vertex vvv where (api_center == 1) do {
					new_spoke := new_edge(vv.id, vvv.id);
					set edge[new_spoke] cbandwidth (band_max + 2);
					set edge[new_spoke] api_spoke 1;
					set edge[new_spoke] e_side bb;
				};
			};

			foreach vertex vv where (cringvertex == (ring_max + 1) && v_side == bb) do {
				foreach vertex vvv where (cringvertex == (ring_max + 1) && v_side == (bb+1)) do {
					new_apiring := new_edge(vv.id, vvv.id);
					set edge[new_apiring] cringedge (ring_max + 1);
					set edge[new_apiring] e_side bb;
				};
			};

			foreach vertex vv where (cringvertex == (ring_max + 1) && v_side == (bb+1)) do {
				foreach vertex vvv where (cringvertex == ring_max && v_side == bb) do {
					new_apidiag := new_edge(vv.id, vvv.id);
					set edge[new_apidiag] e_diag 1;
					set edge[new_apidiag] e_side bb;
					set edge[new_apidiag] bandsq (band_max + 1);
				};
			};
		}; // end of else statement for edges

		}; // end of building edges

		// Loop to build facets
		for ( bb := 1 ; bb < (discret_max+1) ; bb := bb + 1 ) {

			if (bb == discret_max) then {
			// Right triangle
			foreach edge ee1 where (cbandwidth == (band_max + 1) && e_side == 1) do {
				foreach edge ee2 where (bandsq == (band_max + 1) && e_side == bb) do {
					foreach edge ee3 where (cringedge == ring_max && e_side == bb) do {
						new_right_tri := new_facet(ee1.id, ee2.id, ee3.id);
						set facet[new_right_tri] tri (band_max + 1);
						set facet[new_right_tri] f_side bb;
					};
				};
			};

			// Left triangle
			foreach edge ee1 where (bandsq == (band_max + 1) && e_side == bb) do {
				foreach edge ee2 where (cringedge == (ring_max + 1) && e_side == bb) do {
					foreach edge ee3 where (cbandwidth == (band_max + 1) && e_side == bb) do {
						new_left_tri := new_facet(-ee1.id, -ee2.id, -ee3.id);
						set facet[new_left_tri] tri (band_max + 1);
						set facet[new_left_tri] f_side bb;
					};
				};
			};

			// Spoke triangle
			foreach edge ee1 where (api_spoke == 1 && e_side == 1) do {
				foreach edge ee2 where (api_spoke == 1 && e_side == bb) do {
					foreach edge ee3 where (cringedge == (ring_max + 1) && e_side == bb) do {
						new_spoke_tri := new_facet(ee1.id, -ee2.id, ee3.id);
						set facet[new_spoke_tri] tri (band_max + 2);
						set facet[new_spoke_tri] apispoke_tri 1;
						set facet[new_spoke_tri] f_side bb;
					};
				};
			};
			}
			else {
			// Right triangle
			foreach edge ee1 where (cbandwidth == (band_max + 1) && e_side == (bb+1)) do {
				foreach edge ee2 where (bandsq == (band_max + 1) && e_side == bb) do {
					foreach edge ee3 where (cringedge == ring_max && e_side == bb) do {
						new_right_tri := new_facet(ee1.id, ee2.id, ee3.id);
						set facet[new_right_tri] tri (band_max + 1);
						set facet[new_right_tri] f_side bb;
					};
				};
			};

			// Left triangle
			foreach edge ee1 where (bandsq == (band_max + 1) && e_side == bb) do {
				foreach edge ee2 where (cringedge == (ring_max + 1) && e_side == bb) do {
					foreach edge ee3 where (cbandwidth == (band_max + 1) && e_side == bb) do {
						new_left_tri := new_facet(-ee1.id, -ee2.id, -ee3.id);
						set facet[new_left_tri] tri (band_max + 1);
						set facet[new_left_tri] f_side bb;
					};
				};
			};

			// Spoke triangle
			foreach edge ee1 where (api_spoke == 1 && e_side == (bb+1)) do {
				foreach edge ee2 where (api_spoke == 1 && e_side == bb) do {
					foreach edge ee3 where (cringedge == (ring_max + 1) && e_side == bb) do {
						new_spoke_tri := new_facet(ee1.id, -ee2.id, ee3.id);
						set facet[new_spoke_tri] tri (band_max + 2);
						set facet[new_spoke_tri] apispoke_tri 1;
						set facet[new_spoke_tri] f_side bb;
					};
				};
			};
			}; // end of else statement for facets
		}; // end of building facets

	}; // end of 'buapi' for-loop

	band_max := max(edge elsa where (elsa.cbandwidth > 0), elsa.cbandwidth);
	ring_max := max(edge elsa where (elsa.cringedge > 0), elsa.cringedge);
	tri_max := max(facet ff where (ff.tri > 0), ff.tri);

}; // end of 'buildapi' command
buildapi;
dissolve facets where (apispoke_tri == 1);
dissolve edges where (api_spoke == 1);

//****************************************************************************************************
// 'buildpolarapi' only constructs the polar apical cap
buildpolarapi := {

	band_max := max(edge elsa where (elsa.cbandwidth > 0), elsa.cbandwidth);
	ring_max := max(edge elsa where (elsa.cringedge > 0), elsa.cringedge);
	tri_max := max(facet ff where (ff.tri > 0), ff.tri);

	// loop to build vertices, edges, facets on apical pole. Index 'buapi' increases
	// from center to outward periphery direction
	for ( buapi := 1 ; buapi <= api_bnumber ; buapi := buapi + 1 ) {

		r_major := (buapi/api_bnumber)*vscale*apimajor;
		r_minor := (buapi/api_bnumber)*vscale*apiminor;
		angle := 0;

		//***************
		// Build Vertices
		//***************
		// Loop to build vertices at apical pole
		for ( bb := 1 ; bb <= 6*buapi ; bb := bb+1 ) {

			if (buapi == api_bnumber) then {
			angle := (bb * 1/(buapi*6))*2*pi;
			new_apiringvertex := new_vertex(r_major*cos(angle),
	  	                             		r_minor*sin(angle),
											z_height );
			set vertex[new_apiringvertex] v_side bb;
			set vertex[new_apiringvertex] apiringvertex buapi;
			set vertex[new_apiringvertex] cringvertex (ring_max+1);	/* this is the feature we want */
			set vertex[new_apiringvertex] api_ringvangle angle;
			} else {
			angle := (bb * 1/(buapi*6))*2*pi;
			new_apiringvertex := new_vertex(r_major*cos(angle),
	  	                             		r_minor*sin(angle),
											z_height );
			set vertex[new_apiringvertex] v_side bb;
			set vertex[new_apiringvertex] apiringvertex buapi;
			set vertex[new_apiringvertex] api_ringvangle angle;
			};
		};

		// To set constraint on outermost apical polar vertices, to prevent them from moving downwards
		/*
		foreach vertex vv where (apiringvertex == api_bnumber) do {
			set vertex[vv.id] constraint 3;
		};	*/

		//***************
		// Build Edges
		//***************
		// Loop to build apiringedge on apical pole
		for ( bb := 1 ; bb <= 6*buapi ; bb := bb+1 ) {

			if (buapi == api_bnumber) then {
			if (bb == 6*buapi) then {
				// new api ringedge
				foreach vertex vv where (apiringvertex == buapi && v_side == bb) do {
					foreach vertex vvv where (apiringvertex == buapi && v_side == 1) do {
						new_horizontal := new_edge(vv.id, vvv.id);
						set edge[new_horizontal] apiringedge buapi;
						set edge[new_horizontal] cringedge (ring_max+1);
						set edge[new_horizontal] e_side bb;
					};
				};
			}
			else {
				// new apiringedge
				foreach vertex vv where (apiringvertex == buapi && v_side == bb) do {
					foreach vertex vvv where (apiringvertex == buapi && v_side == bb+1) do {
						new_horizontal := new_edge(vv.id, vvv.id);
						set edge[new_horizontal] apiringedge buapi;
						set edge[new_horizontal] cringedge (ring_max+1);
						set edge[new_horizontal] e_side bb;
					};
				};
			};
			}
			else {
			if (bb == 6*buapi) then {
				// new api ringedge
				foreach vertex vv where (apiringvertex == buapi && v_side == bb) do {
					foreach vertex vvv where (apiringvertex == buapi && v_side == 1) do {
						new_horizontal := new_edge(vv.id, vvv.id);
						set edge[new_horizontal] apiringedge buapi;
						set edge[new_horizontal] e_side bb;
					};
				};
			}
			else {
				// new apiringedge
				foreach vertex vv where (apiringvertex == buapi && v_side == bb) do {
					foreach vertex vvv where (apiringvertex == buapi && v_side == bb+1) do {
						new_horizontal := new_edge(vv.id, vvv.id);
						set edge[new_horizontal] apiringedge buapi;
						set edge[new_horizontal] e_side bb;
					};
				};
			}; // end of else for building 'apiringedge'
			};

		}; // end of bb-loop


		// Build edges (vertical, internal slant ones) & facets on apical pole
		if (buapi > 1) then {
			// new api width edges between apical rings that lie along the spoke angles
			foreach vertex vv where (apiringvertex == buapi && (api_ringvangle == 2*pi || api_ringvangle == pi/3 || api_ringvangle == 2*pi/3 || api_ringvangle == pi || api_ringvangle == 4*pi/3 || api_ringvangle == 5*pi/3) ) do {
				set vertex[vv.id] onspoke_vertex 1;

				foreach vertex vvv where (apiringvertex == buapi-1 && api_ringvangle == vv.api_ringvangle ) do {
					new_vertical := new_edge(vvv.id, vv.id);
					set edge[new_vertical] e_side vv.v_side;
					set edge[new_vertical] api_vert 1;
					set edge[new_vertical] apibandwidth (buapi-1);
				};
			};

			// new api width slanted edges between apical rings that are NOT ON the main spoke angles (where 'onspoke_vertex == 0')
			foreach vertex vv where (apiringvertex == buapi && onspoke_vertex == 0) do {
				if (vv.v_side == 1) then {

				foreach vertex vvv where (apiringvertex == buapi-1 && v_side == (buapi-1)*6 ) do {
					new_apislant := new_edge(vvv.id, vv.id);
					set edge[new_apislant] e_side vv.v_side;
					set edge[new_apislant] e_slant 1;
					set edge[new_apislant] apibandsq (buapi-1);
				};
				}
				else {
				vert_num := floor(vv.v_side*(buapi-1)/buapi);

				foreach vertex vvv where (apiringvertex == buapi-1 && v_side == vert_num ) do {
					new_apislant := new_edge(vvv.id, vv.id);
					set edge[new_apislant] e_side vv.v_side;
					set edge[new_apislant] e_slant 1;
					set edge[new_apislant] apibandsq (buapi-1);
				};
				foreach vertex vvv where (apiringvertex == buapi-1 && v_side == vert_num+1 ) do {
					new_apislant := new_edge(vvv.id, vv.id);
					set edge[new_apislant] e_side vv.v_side;
					set edge[new_apislant] e_slant 2;
					set edge[new_apislant] apibandsq (buapi-1);
				};
				};
			};

			foreach vertex vv where (apiringvertex == buapi && v_side == 1) do {
				foreach vertex vvv where (apiringvertex == buapi-1 && v_side == 1 ) do {
					new_apislant := new_edge(vvv.id, vv.id);
					set edge[new_apislant] e_side vv.v_side;
					set edge[new_apislant] e_slant 2;
					set edge[new_apislant] apibandsq (buapi-1);
				};
			};

			//***************
			// Build Facets
			//***************
			// Build internal facets within each apical band
			// the 1st type of facets associated with the api_vert within each apical band
			foreach edge ee1 where (api_vert == 1 && apibandwidth == (buapi-1)) do {
				if (ee1.e_side == buapi*6) then {
				foreach edge ee2 where (apiringedge == buapi && e_side == ee1.e_side) do {
					foreach edge ee3 where (e_slant == 1 && apibandsq == (buapi-1) && e_side == 1) do {
						new_apibandtri := new_facet(ee1.id, ee2.id, -ee3.id);
						set facet[new_apibandtri] apibandtri (buapi-1);
						set facet[new_apibandtri] color yellow;
					};
				};
				}
				else {
				foreach edge ee2 where (apiringedge == buapi && e_side == ee1.e_side) do {
					foreach edge ee3 where (e_slant == 1 && apibandsq == (buapi-1) && e_side == ee2.e_side+1) do {
						new_apibandtri := new_facet(ee1.id, ee2.id, -ee3.id);
						set facet[new_apibandtri] apibandtri (buapi-1);
						set facet[new_apibandtri] color yellow;
					};
				};
				};
			};

			// Build the 2nd type of internal facets within each apical band
			foreach edge ee1 where (e_slant == 1 && apibandsq == (buapi-1)) do {

				if (ee1.e_side == 1) then {
				foreach edge ee2 where (e_slant == 2 && apibandsq == (buapi-1) && e_side == ee1.e_side) do {
					foreach edge ee3 where (apiringedge == (buapi-1) && e_side == (buapi-1)*6 ) do {
						new_apibandtri := new_facet(ee1.id, -ee2.id, -ee3.id);
						set facet[new_apibandtri] apibandtri (buapi-1);
						set facet[new_apibandtri] color yellow;
					};
				};
				}
				else {
				foreach edge ee2 where (e_slant == 2 && apibandsq == (buapi-1) && e_side == ee1.e_side) do {
					foreach edge ee3 where (apiringedge == (buapi-1) && e_side == floor(ee1.e_side*(buapi-1)/buapi) ) do {
						new_apibandtri := new_facet(ee1.id, -ee2.id, -ee3.id);
						set facet[new_apibandtri] apibandtri (buapi-1);
						set facet[new_apibandtri] color yellow;
					};
				};
				};
			};

			// Build the 3rd type of internal facets within each apical band
			foreach edge ee1 where (e_slant == 2 && apibandsq == (buapi-1)) do {
				foreach edge ee2 where (apiringedge == buapi && e_side == ee1.e_side) do {
					foreach edge ee3 where (api_vert == 1 && apibandwidth == (buapi-1) && e_side == ee1.e_side+1) do {
						new_apibandtri := new_facet(ee1.id, ee2.id, -ee3.id);
						set facet[new_apibandtri] apibandtri (buapi-1);
						set facet[new_apibandtri] color yellow;
					};
				};
			};
			// 'highest_vside' is to calculate the maximum value of v_side within each apical ring
			highest_vside := max(vertex vv where (apiringvertex == buapi), vv.v_side);
			foreach edge ee1 where (e_slant == 2 && apibandsq == (buapi-1) && e_side == highest_vside) do {
				foreach edge ee2 where (apiringedge == buapi && e_side == highest_vside) do {
					foreach edge ee3 where (api_vert == 1 && apibandwidth == (buapi-1) && e_side == 1 ) do {
						new_apibandtri := new_facet(ee1.id, ee2.id, -ee3.id);
						set facet[new_apibandtri] apibandtri (buapi-1);
						set facet[new_apibandtri] color yellow;
					};
				};
			};

			// Build the 4th type of internal facets within each apical band
			foreach edge ee1 where (e_slant == 2 && apibandsq == (buapi-1)) do {
				foreach edge ee2 where (apiringedge == buapi && e_side == ee1.e_side) do {
					foreach edge ee3 where (e_slant == 1 && apibandsq == (buapi-1) && e_side == ee1.e_side+1 ) do {
						new_apibandtri := new_facet(ee1.id, ee2.id, -ee3.id);
						set facet[new_apibandtri] apibandtri (buapi-1);
						set facet[new_apibandtri] color yellow;
					};
				};
			};
		}
		else {
		// The following commands to build apical spokes and apical spoke triangles should only be performed for the buapi==1 loop instance
		// which explains why they are positioned in this 'else' statement

			//********************
			// Build Apical Spokes
			//********************
			// build api spoke edges
			foreach vertex vv where (api_center == 1) do {
				foreach vertex vvv where (apiringvertex == buapi) do {
					new_apispoke := new_edge(vv.id, vvv.id);
					set edge[new_apispoke] api_spoke 1;
					set edge[new_apispoke] e_side vvv.v_side;
				};
			};

			//***************
			// Build Facets
			//***************
			// build apical spoke triangles
			// for the 1st to 5th spoke triangles
			foreach edge ee1 where (api_spoke == 1 && e_side != 6) do {
				foreach edge ee2 where (apiringedge == 1 && e_side == ee1.e_side) do {
					foreach edge ee3 where (api_spoke == 1 && e_side == ee1.e_side+1) do {
						new_spoke_tri := new_facet(ee1.id, ee2.id, -ee3.id);
						set facet[new_spoke_tri] apispoke_tri 1;
						set facet[new_spoke_tri] f_side ee1.e_side;
						set facet[new_spoke_tri] color yellow;
					};
				};
			};
			// for the 6th spoke triangle
			foreach edge ee1 where (api_spoke == 1 && e_side == 6) do {
				foreach edge ee2 where (apiringedge == 1 && e_side == ee1.e_side) do {
					foreach edge ee3 where (api_spoke == 1 && e_side == 1) do {
						new_spoke_tri := new_facet(ee1.id, ee2.id, -ee3.id);
						set facet[new_spoke_tri] apispoke_tri 1;
						set facet[new_spoke_tri] f_side ee1.e_side;
						set facet[new_spoke_tri] color yellow;
					};
				};
			};
		};	// end of 'else'

	}; // end of 'buapi' for-loop

	band_max := max(edge elsa where (elsa.cbandwidth > 0), elsa.cbandwidth);
	ring_max := max(edge elsa where (elsa.cringedge > 0), elsa.cringedge);
	tri_max := max(facet ff where (ff.tri > 0), ff.tri);

}; // end of 'buildapi' command
buildpolarapi;

//****************************************************************************************************
// this is only for constructing the new topmost lateral band edges and facets for the subdivision of the topmost lateral band
build_top_lat := {

	// Count/determine status of discretization
	band_max := max(edge elsa where (elsa.cbandwidth > 0), elsa.cbandwidth);
	ring_max := max(edge elsa where (elsa.cringedge > 0), elsa.cringedge);
	tri_max := max(facet ff where (ff.tri > 0), ff.tri);

	// Loop to build edges
	for ( bb := 1 ; bb < (discret_max+1) ; bb := bb + 1 ) {

		if (bb == discret_max) then {
			foreach vertex vv where (cringvertex == (ring_max-1) && v_side == bb) do {
				foreach vertex vvv where (cringvertex == ring_max && v_side == bb) do {
					new_apivertical := new_edge(vv.id, vvv.id);
					set edge[new_apivertical] cbandwidth (band_max + 1);
					set edge[new_apivertical] e_side bb;
				};
			};
			foreach vertex vv where (cringvertex == ring_max && v_side == 1) do {
				foreach vertex vvv where (cringvertex == (ring_max-1) && v_side == bb) do {
					new_apidiag := new_edge(vv.id, vvv.id);
					set edge[new_apidiag] e_diag 1;
					set edge[new_apidiag] e_side bb;
					set edge[new_apidiag] bandsq (band_max + 1);
				};
			};
		}
		else {
			foreach vertex vv where (cringvertex == (ring_max-1) && v_side == bb) do {
				foreach vertex vvv where (cringvertex == ring_max && v_side == bb) do {
					new_apivertical := new_edge(vv.id, vvv.id);
					set edge[new_apivertical] cbandwidth (band_max + 1);
					set edge[new_apivertical] e_side bb;
				};
			};
			foreach vertex vv where (cringvertex == ring_max && v_side == (bb+1)) do {
				foreach vertex vvv where (cringvertex == (ring_max-1) && v_side == bb) do {
					new_apidiag := new_edge(vv.id, vvv.id);
					set edge[new_apidiag] e_diag 1;
					set edge[new_apidiag] e_side bb;
					set edge[new_apidiag] bandsq (band_max + 1);
				};
			};
		}; // end of else statement for edges

		}; // end of building edges

	// Loop to build facets
	for ( bb := 1 ; bb < (discret_max+1) ; bb := bb + 1 ) {

			if (bb == discret_max) then {
			// Right triangle
			foreach edge ee1 where (cbandwidth == (band_max + 1) && e_side == 1) do {
				foreach edge ee2 where (bandsq == (band_max + 1) && e_side == bb) do {
					foreach edge ee3 where (cringedge == (ring_max-1) && e_side == bb) do {
						new_right_tri := new_facet(ee1.id, ee2.id, ee3.id);
						set facet[new_right_tri] tri (band_max + 1);
						set facet[new_right_tri] f_side bb;
					};
				};
			};

			// Left triangle
			foreach edge ee1 where (bandsq == (band_max + 1) && e_side == bb) do {
				foreach edge ee2 where (cringedge == ring_max && e_side == bb) do {
					foreach edge ee3 where (cbandwidth == (band_max + 1) && e_side == bb) do {
						new_left_tri := new_facet(-ee1.id, -ee2.id, -ee3.id);
						set facet[new_left_tri] tri (band_max + 1);
						set facet[new_left_tri] f_side bb;
					};
				};
			};
			}
			else {
			// Right triangle
			foreach edge ee1 where (cbandwidth == (band_max + 1) && e_side == (bb+1)) do {
				foreach edge ee2 where (bandsq == (band_max + 1) && e_side == bb) do {
					foreach edge ee3 where (cringedge == (ring_max-1) && e_side == bb) do {
						new_right_tri := new_facet(ee1.id, ee2.id, ee3.id);
						set facet[new_right_tri] tri (band_max + 1);
						set facet[new_right_tri] f_side bb;
					};
				};
			};

			// Left triangle
			foreach edge ee1 where (bandsq == (band_max + 1) && e_side == bb) do {
				foreach edge ee2 where (cringedge == ring_max && e_side == bb) do {
					foreach edge ee3 where (cbandwidth == (band_max + 1) && e_side == bb) do {
						new_left_tri := new_facet(-ee1.id, -ee2.id, -ee3.id);
						set facet[new_left_tri] tri (band_max + 1);
						set facet[new_left_tri] f_side bb;
					};
				};
			};
			}; // end of else statement for facets
		}; // end of building facets

};
build_top_lat;

//****************************************************************************************************
// Count/determine status of discretization
band_max := max(edge elsa where (elsa.cbandwidth > 0), elsa.cbandwidth);
ring_max := max(edge elsa where (elsa.cringedge > 0), elsa.cringedge);
tri_max := max(facet ff where (ff.tri > 0), ff.tri);

// Do not place the definition of new_body inside a function, it should only be executed/introduced once!
// For example: 'cyl_body := new_body' is what I mean by a new_body definition.
// to set the new facets onto a new center cell body
cyl_body := new_body;

setrefinedbody := {

	set facet frontbody cyl_body where (basespoke_tri == 1);
	set facet frontbody cyl_body where (tri > 0);
	set facet frontbody cyl_body where (apispoke_tri == 1);
	set facet frontbody cyl_body where (apibandtri > 0);

	// Not needed for all vertices of center cell
	//foreach vertex vv where (vv.cringvertex < (ring_max/2)) do {
	//	set vertex[vv.id] constraint 2;
	//};

	// 'frontbody' is used because the orientation of the facets for the neighbouring bands are pointing to the center cell
	set facet frontbody neigh_body where (ntri > 0);
	set facet frontbody neigh_body where (nbase_tri == 1);
	set facet frontbody neigh_body where (napi_tri == 1);
	set facet frontbody neigh_body where (nouter_tri == 1);

	// The purpose of calculating these max index values is to ensure they are updated after every kind of adaptive refinement is executed
	band_max := max(edge elsa where (elsa.cbandwidth > 0), elsa.cbandwidth);
	ring_max := max(edge elsa where (elsa.cringedge > 0), elsa.cringedge);
	tri_max := max(facet ff where (ff.tri > 0), ff.tri);

	nband_max := max(edge elsa where (elsa.nbandwidth > 0), elsa.nbandwidth);
	nring_max := max(edge elsa where (elsa.nringedge > 0), elsa.nringedge);
	ntri_max := max(facet ff where (ff.ntri > 0), ff.ntri);
};
setrefinedbody;
// Remember to call 'setrefinedbody' after every action of refinement (subdivision or coarsening)


// Settings for
// - setting up uniform cortex properties (thickness, depolymerization rate, eqmthickness...)
// the initial value of the signal overactivity time function
sig_t := 1;

// To set initial cortex prop for center cell
// Initial cortex property before basal ring contraction
ini_cortex_prop := {
	set_vol_c := { body[1].target := cvol; };
	set_vol_c;

	band_max := max(edge elsa where (elsa.cbandwidth > 0), elsa.cbandwidth);
	ring_max := max(edge elsa where (elsa.cringedge > 0), elsa.cringedge);
	tri_max := max(facet ff where (ff.tri > 0), ff.tri);

//****************
// For center cell
//****************

	// find the current position of the reference ring which marks the OZ boundary that is a line on a ring
	// Set 'ref_position' with value 0 to indicate OZ begins from basal pole
	ref_position := 0;

	// for base pole facets
	foreach edge eaxbp where (base_spoke == 1 && e_side == 1) do {

		set edge[eaxbp.id] bandthickness ini_bthickness;
		set edge[eaxbp.id] banddepoly ini_bdepoly;
		set edge[eaxbp.id] bandeqmthickness ini_beqmthickness;
		set edge[eaxbp.id] bandviscous ini_bviscous;

		// In this section, the basal pole is not part of OZ, so its zetab is at baseline.
		// zetab is always at baseline level when it is outside of overactive zone signal width
		set edge[eaxbp.id] zetab (zetab_band_baseline);
		eaxbp_tens := (eaxbp.zetab*eaxbp.bandthickness);
		foreach facet eaxialbp where (basespoke_tri == 1 && f_side > 0) do {
			set facet[eaxialbp.id] tension eaxbp_tens;
		};
	};

	// for apical polar cap facets
	foreach edge eaxap where (api_spoke == 1 && e_side == 1) do {

		set edge[eaxap.id] bandthickness ini_bthickness;
		set edge[eaxap.id] banddepoly ini_bdepoly;
		set edge[eaxap.id] bandeqmthickness ini_beqmthickness;
		set edge[eaxap.id] bandviscous ini_bviscous;

		set edge[eaxap.id] ringthickness ini_rthickness;
		set edge[eaxap.id] ringdepoly ini_rdepoly;
		set edge[eaxap.id] ringeqmthickness ini_reqmthickness;
		set edge[eaxap.id] ringviscous ini_rviscous;

		// zetab is at baseline level when it is outside of overactive zone, 2nd degree & higher order away from OZ boundary
		set edge[eaxap.id] zetab (zetab_band_baseline);
		eaxap_tens := (eaxap.zetab*eaxap.bandthickness);
		foreach facet eaxialap where (apispoke_tri == 1 || apibandtri > 0) do {
			set facet[eaxialap.id] tension eaxap_tens;
		};
		set edge[eaxap.id] zetar zetar_ring_baseline;
		eazimu_tens := (eaxap.zetar*eaxap.ringthickness);
		foreach edge eazimu where (cringedge == ring_max && e_side > 0) do {
			set edge[eazimu.id] tension eazimu_tens;
		};
	};

	// for the bands between the 2 poles, not inclusive of the poles
	// If this is a lateral threadmill-style setup, 'band_max' is not at the apical pole. This is more true for the lateral bands of neighbours.
	// As for the center cell, its 'band_max' is indeed at the apical pole, so the 'jojo' loop does not include 'band_max' because the apical pole is already being treated above.
	for ( jojo := 1 ; jojo <= (band_max) ; jojo := jojo + 1 ) {

		// for bands below 'ref_position', such bands are not within the OZ for the extrusion case
		foreach edge eax where (cbandwidth == jojo && cbandwidth < (ref_position) && e_side == 1) do {

			band_position := jojo;

			set edge[eax.id] bandthickness ini_bthickness;
			set edge[eax.id] banddepoly ini_bdepoly;
			set edge[eax.id] bandeqmthickness ini_beqmthickness;
			set edge[eax.id] bandviscous ini_bviscous;

			set edge[eax.id] ringthickness ini_rthickness;
			set edge[eax.id] ringdepoly ini_rdepoly;
			set edge[eax.id] ringeqmthickness ini_reqmthickness;
			set edge[eax.id] ringviscous ini_rviscous;

			// Since the bands below 'ref_position' are always outside of OZ
			set edge[eax.id] zetab (zetab_band_baseline);
			eax_tens := (eax.zetab*eax.bandthickness);
			foreach facet eaxial where (tri == band_position && f_side > 0) do {
				set facet[eaxial.id] tension eax_tens;
			};
			set edge[eax.id] zetar zetar_ring_baseline;
			eazimu_tens := (eax.zetar*eax.ringthickness);
			foreach edge eazimu where (cringedge == band_position && e_side > 0) do {
				set edge[eazimu.id] tension eazimu_tens;
			};
		};

		// for bands at the same level and or above 'ref_position'
		// As long 'jojo' begins from '1', writing 'cbandwidth>=(ref_position) would not mess things up
		foreach edge eax where (cbandwidth == jojo && cbandwidth >= (ref_position) && e_side == 1) do {

			band_position := jojo;

			set edge[eax.id] bandthickness ini_bthickness;
			set edge[eax.id] banddepoly ini_bdepoly;
			set edge[eax.id] bandeqmthickness ini_beqmthickness;
			set edge[eax.id] bandviscous ini_bviscous;

			set edge[eax.id] ringthickness ini_rthickness;
			set edge[eax.id] ringdepoly ini_rdepoly;
			set edge[eax.id] ringeqmthickness ini_reqmthickness;
			set edge[eax.id] ringviscous ini_rviscous;

			// zetab is at baseline level when it is outside of overactive zone, 2nd degree & higher order away from OZ boundary
			set edge[eax.id] zetab (zetab_band_baseline);
			eax_tens := (eax.zetab*eax.bandthickness);
			foreach facet eaxial where (tri == band_position && f_side > 0) do {
				set facet[eaxial.id] tension eax_tens;
			};
			set edge[eax.id] zetar zetar_ring_baseline;
			eazimu_tens := (eax.zetar*eax.ringthickness);
			foreach edge eazimu where (cringedge == band_position && e_side > 0) do {
				set edge[eazimu.id] tension eazimu_tens;
			};
		};

	};

	// for the special magenta basal contracting ring on centre cell
	// On centre cell, to determine the cringedge value of the basal magenta ring
	basal_ring_posit := avg(edge elle where (elle.color == magenta), elle.cringedge);
	// Set attributes and tension on basal contracting ring
	foreach edge bee where (bee.cbandwidth == basal_ring_posit && e_side == 1) do {

		set edge[bee.id] ringthickness basecontract_rthickness;
		set edge[bee.id] ringdepoly basecontract_rdepoly;
		set edge[bee.id] ringeqmthickness basecontract_reqmthickness;
		set edge[bee.id] ringviscous basecontract_rviscous;

		set edge[bee.id] zetar basecontract_zetar_ring;
		eazimu_tens := (bee.zetar*bee.ringthickness);
		foreach edge eazimu where (color == magenta && cringedge == basal_ring_posit && e_side > 0) do {
			set edge[eazimu.id] tension eazimu_tens;
		};
	};

	// coloring	of overactive zone
	foreach edge exx where (e_side == 1 && zetab > ((4/5)*sig_amp) && cbandwidth > 0 ) do {
		band_axx:= exx.cbandwidth;
		foreach facet ff where (ff.tri == band_axx) do {
			set facet[ff.id] color 4;	//red
		};
	};
	foreach edge exx where (e_side == 1 && zetab > ((4/5)*sig_amp) && base_spoke == 1) do {
		foreach facet ff where (ff.basespoke_tri == 1) do {
			set facet[ff.id] color 4;	//red
		};
	};

	foreach edge exx where (e_side == 1 && zetab > ((3/5)*sig_amp) && zetab <= ((4/5)*sig_amp) && cbandwidth > 0 ) do {
		band_axx:= exx.cbandwidth;
		foreach facet ff where (ff.tri == band_axx) do {
			set facet[ff.id] color 12;	//lightred
		};
	};
	foreach edge exx where (e_side == 1 && zetab > ((3/5)*sig_amp) && zetab <= ((4/5)*sig_amp) && base_spoke == 1) do {
		foreach facet ff where (ff.basespoke_tri == 1) do {
			set facet[ff.id] color 12;	//lightred
		};
	};

	foreach edge exx where (e_side == 1 && zetab > ((2/5)*sig_amp) && zetab <= ((3/5)*sig_amp) && cbandwidth > 0 ) do {
		band_axx:= exx.cbandwidth;
		foreach facet ff where (ff.tri == band_axx) do {
			set facet[ff.id] color 2;	//green
		};
	};
	foreach edge exx where (e_side == 1 && zetab > ((2/5)*sig_amp) && zetab <= ((3/5)*sig_amp) && base_spoke == 1) do {
		foreach facet ff where (ff.basespoke_tri == 1) do {
			set facet[ff.id] color 2;	//green
		};
	};

	foreach edge exx where (e_side == 1 && zetab > ((1/5)*sig_amp) && zetab <= ((2/5)*sig_amp) && cbandwidth > 0 ) do {
		band_axx:= exx.cbandwidth;
		foreach facet ff where (ff.tri == band_axx) do {
			set facet[ff.id] color 9;	//lightblue
		};
	};
	foreach edge exx where (e_side == 1 && zetab > ((1/5)*sig_amp) && zetab <= ((2/5)*sig_amp) && base_spoke == 1) do {
		foreach facet ff where (ff.basespoke_tri == 1) do {
			set facet[ff.id] color 9;	//lightblue
		};
	};

	foreach edge exx where (e_side == 1 && zetab > (zetab_band_baseline) && zetab <= ((1/5)*sig_amp) && cbandwidth > 0 ) do {
		band_axx:= exx.cbandwidth;
		foreach facet ff where (ff.tri == band_axx) do {
			set facet[ff.id] color 3;	//cyan
		};
	};
	foreach edge exx where (e_side == 1 && zetab > (zetab_band_baseline) && zetab <= ((1/5)*sig_amp) && base_spoke == 1) do {
		foreach facet ff where (ff.basespoke_tri == 1) do {
			set facet[ff.id] color 3;	//cyan
		};
	};

	foreach edge exx where (e_side == 1 && zetab == (zetab_band_baseline) && cbandwidth > 0 ) do {
		band_axx:= exx.cbandwidth;
		foreach facet ff where (ff.tri == band_axx) do {
			set facet[ff.id] color 14;	// yellow
		};
	};
	foreach edge exx where (e_side == 1 && zetab == (zetab_band_baseline) && base_spoke == 1) do {
		foreach facet ff where (ff.basespoke_tri == 1) do {
			set facet[ff.id] color 14;	// yellow
		};
	};
};
ini_cortex_prop;

//****************************************************************************************************

