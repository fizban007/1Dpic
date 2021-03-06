#include "algorithms/current_deposit_Esirkepov.h"
#include "algorithms/interpolation.h"
#include "utils/util_functions.h"
#include "data/detail/multi_array_utils.hpp"

namespace Aperture {

CurrentDepositer_Esirkepov::CurrentDepositer_Esirkepov(const Environment& env)
    : m_env(env) {}

CurrentDepositer_Esirkepov::~CurrentDepositer_Esirkepov() {}

void CurrentDepositer_Esirkepov::deposit(SimData& data, double dt) {
  Logger::print_detail("Depositing current");
  auto& part = data.particles;
  data.J.initialize();

  for (Index_t i = 0; i < part.size(); i++) {
    data.Rho[i].initialize();
    data.J_s[i].initialize();
    // data.V[i].initialize();
    split_delta_rho(data.J_s[i], data.Rho[i], part[i], dt);
    // normalize_density(data.Rho[i], data.Rho[i]);
  }

  // Handle periodic boundary by copying over the deposited quantities
  if (m_periodic) {
    auto& mesh = data.J.grid().mesh();
    for (int i = 0; i < mesh.guard[0]; i++) {
      // rho
      for (unsigned int j = 0; j < part.size(); j++) {
        data.Rho[j](i + mesh.reduced_dim(0)) += data.Rho[j](i);
        data.Rho[j](i) = 0.0;
        data.Rho[j](2 * mesh.guard[0] - 1 - i) += data.Rho[j](mesh.dims[0] - 1 - i);
        data.Rho[j](mesh.dims[0] - 1 - i) = 0.0;
      }
    }
  }

  // detail::map_multi_array(data.J.data(0), data.J_s[0], data.J.grid().extent(), detail::Op_PlusAssign<Scalar>());
  // communication on the just deposited Rho
  if (m_comm_rho != nullptr) {
    for (Index_t i = 0; i < part.size(); i++) {
      m_comm_rho(data.Rho[i]);
    }
  }
  // Now we have delta Q in every cell, add them up along all directions

  scan_current(data.J_s[0]);
  scan_current(data.J_s[1]);
  detail::map_multi_array(data.J.data(0).begin(), data.J_s[0].data().begin(),
                          data.J.grid().extent(), detail::Op_PlusAssign<Scalar>());
  detail::map_multi_array(data.J.data(0).begin(), data.J_s[1].data().begin(),
                          data.J.grid().extent(), detail::Op_PlusAssign<Scalar>());
  // for (unsigned int j = 0; j < part.size(); j++) {
  //   normalize_velocity(data.Rho[j], data.V[j]);
  // }
  // Call communication on just scanned J
  if (m_comm_J != nullptr) {
    m_comm_J(data.J);
  }

  // periodic boundary issues
  auto& mesh = data.J.grid().mesh();
  if (m_periodic) {
    for (int i = 0; i < mesh.guard[0]; i++) {
      data.J(0, i + mesh.reduced_dim(0)) += data.J(0, i);
      data.J(0, i) = 0.0;
      data.J(0, 2 * mesh.guard[0] - 1 - i) += data.J(0, mesh.dims[0] - 1 - i);
      data.J(0, mesh.dims[0] - 1 - i) = 0.0;
    }
    data.J(0, mesh.guard[0] - 1) = data.J(0, mesh.reduced_dim(0) + mesh.guard[0] - 1);
  }
}

void CurrentDepositer_Esirkepov::split_delta_rho(vfield& J, sfield& Rho,
                                                 const Particles& particles,
                                                 double dt) {
  Interpolator interp(m_interp);
  auto& part = particles.data();
  auto& grid = J.grid();
  auto charge = particles.charge();
  if (grid.dim() == 1) {
    // Logger::print_info("Computing rho");
    // loop over all the particles
    for (Index_t n = 0; n < particles.number(); n++) {
      if (particles.is_empty(n)) continue;
      // double v = part.p1[n] / part.gamma[n];

      int c = part.cell[n];
      int c_p = c;
      auto x = part.x1[n];
      auto x_p = part.x1[n] - part.dx1[n];
      c_p += std::floor(x_p);
      x_p -= (double)c_p - c;
      // Logger::print_info("{}, {}, {}, {}", c, c_p, x, x_p);
      // auto v3 = part.dx3[n];

      double s0, s1;
      for (int i = c_p - interp.radius() - 1;
           i <= c_p + interp.support() - interp.radius(); i++) {
        double w0 = interp.interp_cell(x, c, i);
        s1 = w0;
        // This is kinda redundant
        // int idx = grid.mesh().get_idx(i);
        // assert(idx < grid.mesh().dims[0] && idx >= 0);

        if (!check_bit(part.flag[n], ParticleFlag::ignore_current)) {
          double w0_p = interp.interp_cell(x_p, c_p, i);
          // double w1_p = interp.interp_cell(x_p[1], c_p[1], j);
          s0 = w0_p;
          // double s00 = w0_p * w1_p;
          J(0, i) += -charge * (s1 - s0) * grid.mesh().delta[0] / dt;
        }
        Rho(i) += charge * s1;
        // V(i) += v * s1;
        // Logger::print_info("weights are {}, {}; {}", s0, s1, s1-s0);
      }
    }
  }
}

void CurrentDepositer_Esirkepov::split_delta_rho(sfield& J, sfield& Rho,
                                                 const Particles& particles,
                                                 double dt) {
  Interpolator interp(m_interp);
  auto& part = particles.data();
  auto& grid = J.grid();
  auto charge = particles.charge();
  if (grid.dim() == 1) {
    // Logger::print_info("Computing rho");
    // loop over all the particles
    for (Index_t n = 0; n < particles.number(); n++) {
      if (particles.is_empty(n)) continue;
      // double v = part.p1[n] / part.gamma[n];

      int c = part.cell[n];
      int c_p = c;
      auto x = part.x1[n];
      auto x_p = part.x1[n] - part.dx1[n];
      c_p += std::floor(x_p);
      x_p -= (double)c_p - c;
      // Logger::print_info("{}, {}, {}, {}", c, c_p, x, x_p);
      // auto v3 = part.dx3[n];

      double s0, s1;
      for (int i = c_p - interp.radius() - 1;
           i <= c_p + interp.support() - interp.radius(); i++) {
        double w0 = interp.interp_cell(x, c, i);
        s1 = w0;
        // This is kinda redundant
        // int idx = grid.mesh().get_idx(i);
        // assert(idx < grid.mesh().dims[0] && idx >= 0);

        if (!check_bit(part.flag[n], ParticleFlag::ignore_current)) {
          double w0_p = interp.interp_cell(x_p, c_p, i);
          // double w1_p = interp.interp_cell(x_p[1], c_p[1], j);
          s0 = w0_p;
          // double s00 = w0_p * w1_p;
          J(i) += -charge * (s1 - s0) * grid.mesh().delta[0] / dt;
        }
        Rho(i) += charge * s1;
        // Logger::print_info("weights are {}, {}; {}", s0, s1, s1-s0);
      }
    }
  }
}

void CurrentDepositer_Esirkepov::scan_current(sfield& J) {
  auto& grid = J.grid();
  if (grid.dim() == 1) {
    for (int i = 1; i < grid.mesh().dims[0]; i++) {
      J.data()[i] += J.data()[i-1];
    }
  }
}

void CurrentDepositer_Esirkepov::scan_current(vfield& J) {
  auto& grid = J.grid();
  if (grid.dim() == 1) {
    for (int i = 1; i < grid.mesh().dims[0]; i++) {
      J.data(0)[i] += J.data(0)[i-1];
    }
  }
}

// FIXME: Boundary conditions!
void CurrentDepositer_Esirkepov::normalize_current(const vfield& I, vfield& J) {
  auto& grid = I.grid();
  auto& mesh = I.grid().mesh();
  if (grid.dim() == 1) {
    for (int i = 0; i < mesh.dims[0]; i++) {
      J.data(0)[i] = I.data(0)[i] * mesh.delta[0];
    }
  }
}

void CurrentDepositer_Esirkepov::normalize_density(const sfield& Q,
                                                   sfield& rho, sfield& V) {
}

void CurrentDepositer_Esirkepov::normalize_velocity(const sfield& rho, sfield& V) {
  auto& grid = rho.grid();
  auto& mesh = grid.mesh();
  if (grid.dim() == 1) {
    for (int i = 0; i < mesh.dims[0]; i++) {
      // J.data(0)[i] = I.data(0)[i] * mesh.delta[0];
      if (std::abs(rho(i)) > 1e-5)
        V(i) /= rho(i);
      else
        V(i) = 0.0;
    }
  }
}

}
