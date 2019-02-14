const N = 32,
	SCALE = 10;

var t = 0,
	iter = 1;

function setup() {
	cnv = createCanvas(N * SCALE, N * SCALE, P2D);
	fluid = new FluidCube(1, 0.2, 0, 0.0000001);

	// label = createDiv('Iterations');
	// sliderIt = createSlider(0, 64, 1);
	// sliderIt.parent(label);
}

function draw() {
	background(0);
	var cx = (0.5 * width) / SCALE;
	var cy = (0.5 * height) / SCALE;

	fluid.addDensity(cx, cy, Math.floor(Math.random() * 150) + 50);

	var deltaX = mouseX / SCALE - cx;
	var deltaY = mouseY / SCALE - cy;
	var rad = Math.atan2(deltaY, deltaX);

	console.log(mouseX, mouseY, cx, cy);

	// var angle = rad * (180 / PI);
	var v = p5.Vector.fromAngle(rad);
	v.mult(0.2);
	t += 0.01;
	fluid.addVelocity(cx, cy, v.x, v.y);

	fluid.step();
	fluid.render();
	fluid.fade();

	// iter = sliderIt.value();
}

function IX(x, y) {
	x = constrain(x, 0, N - 1);
	y = constrain(y, 0, N - 1);
	return x + y * N;
}

function diffuse(b, x, x0, diff, dt, iter, N) {
	var a = dt * diff * (N - 2) * (N - 2);
	lin_solve(b, x, x0, a, 1 + 6 * a, iter, N);
}

function lin_solve(b, x, x0, a, c, iter, N) {
	var cRecip = 1.0 / c;
	for (var k = 0; k < iter; k++) {
		for (var j = 1; j < N - 1; j++) {
			for (var i = 1; i < N - 1; i++) {
				x[IX(i, j)] =
					(x0[IX(i, j)] +
						a *
							(x[IX(i + 1, j)] +
								x[IX(i - 1, j)] +
								x[IX(i, j + 1)] +
								x[IX(i, j - 1)])) *
					cRecip;
			}
		}

		set_bnd(b, x, N);
	}
}

function project(velocX, velocY, p, div, iter, N) {
	for (var j = 1; j < N - 1; j++) {
		for (var i = 1; i < N - 1; i++) {
			div[IX(i, j)] =
				(-0.5 *
					(velocX[IX(i + 1, j)] -
						velocX[IX(i - 1, j)] +
						velocY[IX(i, j + 1)] -
						velocY[IX(i, j - 1)])) /
				N;
			p[IX(i, j)] = 0;
		}
	}
	set_bnd(0, div, N);
	set_bnd(0, p, N);
	lin_solve(0, p, div, 1, 6, iter, N);

	for (var j = 1; j < N - 1; j++) {
		for (var i = 1; i < N - 1; i++) {
			velocX[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * N;
			velocY[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * N;
		}
	}

	set_bnd(1, velocX, N);
	set_bnd(2, velocY, N);
}

function advect(b, d, d0, velocX, velocY, dt, N) {
	var i0, i1, j0, j1;

	var dtx = dt * (N - 2);
	var dty = dt * (N - 2);

	var s0, s1, t0, t1;
	var tmp1, tmp2, x, y;

	var Nfloat = N;
	var ifloat, jfloat;
	var i, j;

	for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
		for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
			tmp1 = dtx * velocX[IX(i, j)];
			tmp2 = dty * velocY[IX(i, j)];
			x = ifloat - tmp1;
			y = jfloat - tmp2;

			if (x < 0.5) x = 0.5;
			if (x > Nfloat + 0.5) x = Nfloat + 0.5;
			i0 = Math.floor(x);
			i1 = i0 + 1.0;
			if (y < 0.5) y = 0.5;
			if (y > Nfloat + 0.5) y = Nfloat + 0.5;
			j0 = Math.floor(y);
			j1 = j0 + 1.0;

			s1 = x - i0;
			s0 = 1.0 - s1;
			t1 = y - j0;
			t0 = 1.0 - t1;

			var i0i = i0;
			var i1i = i1;
			var j0i = j0;
			var j1i = j1;

			d[IX(i, j)] =
				s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)]) +
				s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);
		}
	}

	set_bnd(b, d, N);
}

function set_bnd(b, x, N) {
	for (var k = 1; k < N - 1; k++) {
		for (var i = 1; i < N - 1; i++) {
			x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
			x[IX(i, N - 1)] = b == 2 ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
		}
	}
	for (var k = 1; k < N - 1; k++) {
		for (var j = 1; j < N - 1; j++) {
			x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
			x[IX(N - 1, j)] = b == 1 ? -x[IX(N - 2, j)] : x[IX(N - 2, j)];
		}
	}

	x[IX(0, 0)] = 0.33 * (x[IX(1, 0, 0)] + x[IX(0, 1, 0)] + x[IX(0, 0, 1)]);
	x[IX(0, N - 1)] = 0.33 * (x[IX(1, N - 1, 0)] + x[IX(0, N - 2, 0)] + x[IX(0, N - 1, 1)]);
	x[IX(N - 1, 0)] = 0.33 * (x[IX(N - 2, 0, 0)] + x[IX(N - 1, 1, 0)] + x[IX(N - 1, 0, 1)]);
	x[IX(N - 1, N - 1)] =
		0.33 * (x[IX(N - 2, N - 1, 0)] + x[IX(N - 1, N - 2, 0)] + x[IX(N - 1, N - 1, 1)]);
}
