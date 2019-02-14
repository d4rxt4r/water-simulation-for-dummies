class FluidCube {
	constructor(size, dt, diffusion, viscosity) {
		this.size = size;
		this.dt = dt;
		this.diff = diffusion;
		this.visc = viscosity;

		this.s = new Array(N * N).fill(0);
		this.density = new Array(N * N).fill(0);

		this.Vx = new Array(N * N).fill(0);
		this.Vy = new Array(N * N).fill(0);

		this.Vx0 = new Array(N * N).fill(0);
		this.Vy0 = new Array(N * N).fill(0);
	}

	addDensity(x, y, amount) {
		this.density[IX(x, y)] += amount;
	}

	addVelocity(x, y, amountX, amountY) {
		var index = IX(x, y);

		this.Vx[index] += amountX;
		this.Vy[index] += amountY;
	}

	step() {
		var visc = this.visc,
			diff = this.diff,
			dt = this.dt,
			Vx = this.Vx,
			Vy = this.Vy,
			Vx0 = this.Vx0,
			Vy0 = this.Vy0,
			s = this.s,
			density = this.density;

		diffuse(1, Vx0, Vx, visc, dt, iter, N);
		diffuse(2, Vy0, Vy, visc, dt, iter, N);

		project(Vx0, Vy0, Vx, Vy, iter, N);

		advect(1, Vx, Vx0, Vx0, Vy0, dt, N);
		advect(2, Vy, Vy0, Vx0, Vy0, dt, N);

		project(Vx, Vy, Vx0, Vy0, iter, N);

		diffuse(0, s, density, diff, dt, iter, N);
		advect(0, density, s, Vx, Vy, dt, N);
	}

	render() {
		colorMode(HSB, 255);

		for (var i = 0; i < N; i++) {
			for (var j = 0; j < N; j++) {
				var x = i * SCALE;
				var y = j * SCALE;
				var d = this.density[IX(i, j)];
				fill(d, 255);
				// fill(255, 100);
				// noStroke();
				stroke('black');
				rect(x, y, SCALE, SCALE);
			}
		}
  }
  
  fade() {
    for (var i = 0; i < this.density.length; i++) {
      var d = this.density[i];
      this.density[i] = constrain(d - 0.02, 0, 255);
    }
  }
}
