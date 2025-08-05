import argparse
import math

def main():
    parser = argparse.ArgumentParser(description="Generate SBEND-based wiggler for Elegant")
    parser.add_argument("-K", type=float, required=True, help="Wiggler parameter K")
    parser.add_argument("-energy", type=float, required=True, help="Beam energy in MeV")
    parser.add_argument("-period", type=float, required=True, help="Wiggler period in cm")
    parser.add_argument("-length", type=float, required=True, help="Total wiggler length in meters")
    parser.add_argument("-rootname", type=str, required=True, help="Base name for SBEND elements")
    parser.add_argument("-name", type=str, default="WIGGLER", help="Name of the wiggler LINE")
    args = parser.parse_args()

    # Convert units
    period_m = args.period / 100.0               # cm → m
    gamma = args.energy / 0.511                  # Lorentz factor
    beam_energy_GeV = args.energy / 1000.0
    brho = beam_energy_GeV / 0.299792458         # Magnetic rigidity [T·m]

    # Calculate peak magnetic field from K
    B = args.K / (93.4 * period_m)               # Tesla
    print("B = ", B)

    # Approximate number of dipoles: 2 per period
    half_period = period_m / 2.0
    num_dipoles = round(args.length / half_period)
    dipole_length = args.length / num_dipoles    # m
    print("Dipole_length = ", dipole_length)
    angle = B * dipole_length / brho             # radians
    print("Angle = ", angle)

    with open(f"{args.rootname}.lte", "w") as f:
        f.write(f"! Generated wiggler: {args.K=}, {args.energy=} MeV, {args.period=} cm, {args.length=} m\n")
        for i in range(num_dipoles):
            sign = 1 if i % 2 == 0 else -1
            label = f"{args.rootname.upper()}{i:03d}"
            f.write(f"{label}: CSRCSBEND, N_KICKS=40,INTEGRATION_ORDER=4,ISR=1,SYNCH_RAD=1,BINS=1500,CSR = 1, SG HALFWIDTH = 2, L = {dipole_length:.6f}, ANGLE = {sign * angle:.6f}\n")

        f.write(f"\n{args.name}: LINE = (")
        f.write(", ".join(f"{args.rootname.upper()}{i:03d}" for i in range(num_dipoles)))
        f.write(")\n")

    print(f"Wrote wiggler lattice to '{args.rootname}.lte' with {num_dipoles} CSRCSBENDs")

if __name__ == "__main__":
    main()