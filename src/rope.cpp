#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass,
           float k, vector<int> pinned_nodes) {
  // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and
  // containing `num_nodes` nodes.
  //
  // New masses can be created (see mass.h):
  // Mass *m = new Mass(position, mass, bool)
  // and should be added to the masses vector.
  //
  // Springs can be created (see spring.h):
  // Spring *s = new Spring(mass_a, mass_b, spring_const)
  // and should be added to the springs vector.
  //
  // Masses corresponding to indices in pinned_nodes
  // should have their pinned field set to true.

  Vector2D offset = (end-start)/(num_nodes-1);

  for(int i = 0; i < num_nodes; i += 1) {
      Vector2D pos = start + (offset*i);
      Mass* m = new Mass(pos, node_mass, false);
      masses.push_back(m);

      if (i > 0){
          Mass* m0 = masses[i-1];
          Spring* s = new Spring(m0, m, k);
          springs.push_back(s);
      }
  }

  for(int i : pinned_nodes) {
      masses[i]->pinned = true;
  }

}

void Rope::simulateEuler(float delta_t, Vector2D gravity) {
  for (auto &s : springs) {
    // TODO (Part 2.1): Use Hooke's law to calculate the force on a node
    Vector2D dist = s->m2->position - s->m1->position;
    float deformation = dist.norm() - s->rest_length;

    Vector2D f_ab = s->k * dist.unit() * deformation;

    // TODO (Part 4.1): Add damping forces
    Vector2D rel_velocity = s->m2->velocity - s->m1->velocity;
    Vector2D damping = 0.5 * rel_velocity;

    // TODO Apply forces as appropriate.
    s->m1->forces += f_ab + damping;
    s->m2->forces -= f_ab + damping;

  }

  for (auto &m : masses) {
    if (!m->pinned) {
      // TODO (Part 2.1): Add the force due to gravity, then compute the new
      // velocity and position

      // f = ma -> f/m = a
      m->velocity += ((m->forces/m->mass)+gravity)*delta_t;
      m->position += m->velocity*delta_t;

    }

    // TODO Reset all forces on each mass
    m->forces = Vector2D();

  }
}

void Rope::simulateVerlet(float delta_t, Vector2D gravity) {
// TODO (Part 3.1): Clear forces
  for(Mass* m : masses) {
      m->forces = Vector2D();
  }
// TODO (Part 3.1): Simulate one timestep of the rope using explicit Verlet
  for (auto &s : springs) {
    Vector2D dist = s->m2->position - s->m1->position;
    float deform = dist.norm() - s->rest_length;

    Vector2D f_ab = s->k * dist.unit() * deform;
    s->m1->forces += f_ab;
    s->m2->forces -= f_ab;

  }

  for (auto &m : masses) {
    if (!m->pinned) {
      Vector2D temp_position = m->position;

      // TODO (Part 3.1): Set the new position of the rope mass
      float offset = (1-0.0005);

      m->position = temp_position + offset*(temp_position - m->last_position) + (((m->forces/m->mass)+gravity)*delta_t*delta_t);
      m->last_position = temp_position;
      // TODO (Part 4.2): Add global Verlet damping

    }
  }
}
}
