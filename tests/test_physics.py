import numpy as np
import pytest
import rebound
from keeler.encounter import Encounter, integrate, elements, jacobi, HILL, OMEGA
from keeler.grains import Box, setup, state, initial_particles
from keeler.optics import slab, shadow_length


def test_free_shear():
    b = np.array([-25e3, 20e3])
    t, q = integrate(b, Encounter(mass_factor=0, samples=51))
    np.testing.assert_allclose(q[:,:,0], np.broadcast_to(b/HILL, (len(t),2)), atol=2e-9)
    np.testing.assert_allclose(q[:,:,1], -1.5*t[:,None]*b/HILL, atol=2e-9)


def test_jacobi_and_guiding_center_shift():
    t,q = integrate(np.array([15e3,20e3,25e3]), Encounter(samples=501))
    c = jacobi(q)
    assert np.max(abs(c-c[0])) < 1e-7
    guiding, radial, vertical = elements(q)
    assert np.all(guiding[-1] > guiding[0])
    r0 = np.linalg.norm(q[0,:,:3],axis=1)
    rf = np.linalg.norm(q[-1,:,:3],axis=1)
    np.testing.assert_allclose(.75*(guiding[-1]**2-guiding[0]**2),
                               radial[-1]**2+6/r0-6/rf, atol=2e-7)


def test_mirror_symmetry():
    _,q=integrate(np.array([-20e3,20e3]), Encounter(samples=101))
    np.testing.assert_allclose(q[:,0],-q[:,1], atol=1e-10)


def test_coplanar_moon_does_not_raise_planar_test_particles():
    _,q=integrate(np.array([20e3]), Encounter(samples=51))
    assert np.max(abs(q[:,:,2])) == 0
    _,q=integrate(np.array([20e3]), Encounter(vertical_amplitude_m=8600, samples=51))
    assert np.max(abs(q[:,:,2]))*HILL > 10


def test_sei_free_vertical_energy():
    c=Box(size_m=32, radius_m=2, self_gravity=False, collisions=False)
    sim,callback=setup(c)
    q0=state(sim)
    sim.steps(1000)
    q=state(sim)
    np.testing.assert_allclose(q[:,5]**2+OMEGA**2*q[:,2]**2,
                               q0[:,5]**2+OMEGA**2*q0[:,2]**2, rtol=1e-10)


@pytest.mark.parametrize('epsilon',[0.,.5,1.])
def test_binary_collision_momentum_and_energy(epsilon):
    s=rebound.Simulation()
    s.gravity='none'
    s.integrator='leapfrog'
    s.collision='direct'
    s.collision_resolve='hardsphere'
    s.dt=.005
    def cor(pointer,speed): return epsilon
    s.coefficient_of_restitution=cor
    s.add(m=2,r=1,x=-1.1,vx=1)
    s.add(m=2,r=1,x=1.1,vx=-1)
    s.steps(25)
    assert abs(sum(p.m*p.vx for p in s.particles)) < 1e-12
    assert abs(s.particles[0].vx+epsilon) < 1e-12
    assert abs(sum(.5*p.m*p.vx**2 for p in s.particles)-2*epsilon**2) < 1e-12


def test_initial_nonoverlap_and_physical_scaling():
    a,b=Box(radius_m=1),Box(radius_m=2)
    assert a.count/b.count == pytest.approx(8,rel=.02)
    assert a.optical_depth/b.optical_depth == pytest.approx(2,rel=.02)
    q=initial_particles(b)
    for i in range(len(q)):
        d=q[i+1:,:3]-q[i,:3]
        d[:,:2]-=b.size_m*np.round(d[:,:2]/b.size_m)
        assert np.all(np.sum(d*d,axis=1)>=(2*b.radius_m)**2)


def test_optics_limits():
    assert slab(0)==0
    assert slab(.3)<slab(.6)<slab(100)
    assert shadow_length(1000,1)==pytest.approx(57289.96163)
    with pytest.raises(ValueError): slab(-1)
    with pytest.raises(ValueError): shadow_length(10,0)


def test_tree_gravity_close_to_direct():
    base,_=setup(Box(size_m=64, self_gravity=True,collisions=False))
    base.N_ghost_x=base.N_ghost_y=0
    base.gravity='basic'
    base.steps(1)
    ref=np.array([[p.ax,p.ay,p.az] for p in base.particles])
    tree,_=setup(Box(size_m=64,self_gravity=True,collisions=False))
    tree.N_ghost_x=tree.N_ghost_y=0
    tree.opening_angle2=.1**2
    tree.steps(1)
    actual=np.array([[p.ax,p.ay,p.az] for p in tree.particles])
    assert np.linalg.norm(ref)>1e-8
    assert np.linalg.norm(actual-ref)/np.linalg.norm(ref)<.01


def test_vertical_extent_does_not_wrap_splash():
    config=Box(size_m=32, radius_m=2, self_gravity=False, collisions=False)
    sim,callback=setup(config)
    sim.particles[0].z=50
    sim.particles[0].vz=0
    sim.steps(1)
    assert sim.particles[0].z == pytest.approx(50*np.cos(OMEGA*config.dt_s),abs=1e-10)


def test_oblique_collision_redirects_horizontal_energy_vertically():
    s=rebound.Simulation()
    s.gravity='none';s.integrator='leapfrog'
    s.collision='direct';s.collision_resolve='hardsphere';s.dt=.001
    def cor(pointer,speed): return .5
    s.coefficient_of_restitution=cor
    s.add(m=1,r=1,x=-1,z=-.5,vx=1)
    s.add(m=1,r=1,x=1,z=.5,vx=-1)
    s.steps(150)
    assert abs(s.particles[0].vz)>.1
    assert abs(sum(p.m*p.vz for p in s.particles))<1e-12
    assert sum(.5*p.m*(p.vx**2+p.vz**2) for p in s.particles)<1


def test_crowding_of_unforced_shear_is_unity():
    from keeler.encounter import crowding
    b=np.linspace(15e3,25e3,11)
    t,q=integrate(b,Encounter(mass_factor=0,samples=31))
    ratio,det=crowding(t,b,q)
    np.testing.assert_allclose(ratio,1,atol=2e-12)
    assert (det<0).all()


def test_step_bound_rejects_unresolved_grain_run():
    from keeler.grains import run
    with pytest.raises(RuntimeError,match='Step too large'):
        run(Box(dt_s=100),samples=2)


def test_contact_check_includes_sheared_images():
    from keeler.contacts import deep_overlap_fraction
    from keeler.encounter import PERIOD
    q=np.zeros((2,6))
    q[0,:2]=[-63,0]
    q[1,:2]=[63,0]
    assert deep_overlap_fraction(q,128,2,0)==1
    # After a quarter orbit the neighboring radial box is shifted in azimuth.
    assert deep_overlap_fraction(q,128,2,PERIOD/4)==0


def test_bridges_law_uses_meters_per_second():
    from keeler.grains import coefficient
    assert coefficient('bridges',.01)==pytest.approx(.32)
    assert coefficient('bridges',-.01)==pytest.approx(.32)
    assert coefficient('bridges',0)==1
    assert coefficient('bridges',1e-10)==1
    assert coefficient('bridges',.1)<.32
    with pytest.raises(ValueError): Box(restitution='unknown')


def test_box_collision_switch_controls_vertical_heating():
    config=Box(size_m=32,orbits=.2)
    collisional,callback=setup(config)
    free,free_callback=setup(Box(size_m=32,orbits=.2,collisions=False))
    collisional.steps(5000)
    free.steps(5000)
    a,b=state(collisional),state(free)
    ea=np.mean(a[:,5]**2+OMEGA**2*a[:,2]**2)
    eb=np.mean(b[:,5]**2+OMEGA**2*b[:,2]**2)
    assert ea>2*eb
