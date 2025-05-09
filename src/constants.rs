
// Vocoder
pub const FRAMEPERIOD: f32 = 2.0;

pub const CHANNELS_IN_OCTAVE: f32 = 5.0;
pub const SKIP_RATIO: usize = 2;
pub const FLOOR_F0_W: f32 = 80.0;
pub const CEIL_F0_W: f32 = 720.0;
pub const FFT: usize = 1024;
pub const FFT2: usize = 2048;
pub const FFT3: usize = 4096;

// World
pub const FLOOR_F0: f32 = 71.0;
pub const CEIL_F0: f32 = 880.0;
pub const DEFAULT_F0: f32 = 500.0;

// StoneMask
pub const MY_SAFE_GUARD_MINIMUM: f32 = 1e-12; 
pub const FLOOR_F0_STONE_MASK: f32 = 40.0;

// Platinum
pub const LOWEST_F0: f32 = 40.0;

// Matlab
pub const NMIN: usize = 256;

// Comb Sort
pub const SHRINK: f32 = 1.2473309501039787;

// Math
pub const LN_2: f32 = 0.693147180559945309417232121458176568075500134360255254120680009493393621969694715605863326996418687;
pub const LOG2_E: f32 = 1.4426950408889634073599246810018921374266459541530;
pub const PI_HALF: f32 = 1.57079632679489661923132169163975144209858469968755291048747229615390820314310449931401741267105853;
pub const PI: f32 = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214;
pub const TAU: f32 = 6.28318530717958647692528676655900576839433879875021164194988918461563281257241799725606965068423413;
