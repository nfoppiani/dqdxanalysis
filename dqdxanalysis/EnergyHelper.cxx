#ifndef ENERGYHELPER_CXX
#define ENERGYHELPER_CXX

#include "EnergyHelper.h"

namespace lee
{

EnergyHelper::EnergyHelper()
{
  detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
}

void EnergyHelper::energyFromHits(recob::PFParticle const &pfparticle,
                                  std::vector<int> &nHits,
                                  std::vector<double> &pfenergy,
                                  const art::Event &evt,
                                  std::string _pfp_producer)
{
  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);

  art::FindManyP<recob::Cluster> cluster_pfparticle_ass(pfparticle_handle, evt,
                                                        _pfp_producer);
  std::vector<art::Ptr<recob::Cluster>> clusters =
      cluster_pfparticle_ass.at(pfparticle.Self());



  nHits.resize(3);
  pfenergy.resize(3);

  for (size_t icl = 0; icl < clusters.size(); icl++)
  {
    art::FindManyP<recob::Hit> hit_cluster_ass(cluster_handle, evt,
                                               _pfp_producer);
    std::vector<art::Ptr<recob::Hit>> hits =
        hit_cluster_ass.at(clusters[icl].key());

    for (size_t ihit = 0; ihit < hits.size(); ++ihit)
    {
      auto plane_nr = hits[ihit]->View();
      if (plane_nr > 2 || plane_nr < 0)
        continue;

      pfenergy[plane_nr] += hits[ihit]->Integral() * _gain[plane_nr] * work_function / recombination_factor / 1000; // convert MeV to GeV
      nHits[plane_nr]++;
    }
  }
}

double EnergyHelper::trackEnergy_dedx(const art::Ptr<recob::Track> &track,
                                      const art::Event &evt,
                                      std::string _pfp_producer)
{
  auto const &track_handle = evt.getValidHandle<std::vector<recob::Track>>(_pfp_producer);
  art::FindManyP<anab::Calorimetry> calo_track_ass(track_handle, evt, "pandoraNucalo");

  const std::vector<art::Ptr<anab::Calorimetry>> calos = calo_track_ass.at(track->ID());

  double E = 0;

  for (size_t ical = 0; ical < calos.size(); ++ical)
  {
    if (E != 0)
      continue;
    if (!calos[ical])
      continue;

    if (!calos[ical]->PlaneID().isValid)
      continue;

    int planenum = calos[ical]->PlaneID().Plane;

    if (planenum < 0 || planenum > 2)
      continue;
    if (planenum != 2)
      continue;                              // Use informartion from collection plane only
    E = calos[ical]->KineticEnergy() / 1000; // convert to GeV
  }
  return E;
}

void EnergyHelper::nHits(size_t pfp_id,
                         const art::Event &evt,
                         std::vector<int> &nHits,
                         std::string _pfp_producer)
{

  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);

  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt, _pfp_producer);

  std::vector<art::Ptr<recob::Cluster>> clusters = clusters_per_pfpart.at(pfp_id);
  nHits.resize(3);

  art::FindManyP<recob::Hit> hits_per_clusters(cluster_handle, evt, _pfp_producer);

  for (size_t icl = 0; icl < clusters.size(); icl++)
  {
    nHits[clusters[icl]->Plane().Plane] = hits_per_clusters.at(clusters[icl].key()).size();
  }
}

void EnergyHelper::PCA(size_t pfp_id,
                       const art::Event &evt,
                       std::vector<std::vector<double>> &pca_planes,
                       std::string _pfp_producer)
{

  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);

  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt, _pfp_producer);

  std::vector<art::Ptr<recob::Cluster>> clusters = clusters_per_pfpart.at(pfp_id);

  art::FindManyP<recob::Hit> hits_per_clusters(cluster_handle, evt, _pfp_producer);

  double drift = detprop->DriftVelocity() * 1e-3;
  double fromTickToNs = 4.8 / detprop->ReadOutWindowSize() * 1e6;
  double wireSpacing = 0.3;

  for (size_t icl = 0; icl < clusters.size(); icl++)
  {
    TPrincipal fPrincipal(2, "D");

    std::vector<art::Ptr<recob::Hit>> hits = hits_per_clusters.at(clusters[icl].key());

    for (auto &hit : hits)
    {
      double data[2];
      double w = hit->WireID().Wire * wireSpacing;
      double t = fromTickToNs * drift * hit->PeakTime();
      data[0] = w;
      data[1] = t;
      fPrincipal.AddRow(data);
    }

    fPrincipal.MakePrincipals();
    pca_planes[clusters[icl]->Plane().Plane][0] = (*fPrincipal.GetEigenValues())[0];
    pca_planes[clusters[icl]->Plane().Plane][1] = (*fPrincipal.GetEigenValues())[1];
  }
}

void EnergyHelper::dQdx(size_t pfp_id,
                            const art::Event &evt,
                            std::vector<double> &dqdx,
                            std::vector<double> &dqdx_hits,
                            std::vector<int> &dqdx_wires,
                            std::vector<double> &box_start,
                            std::vector<double> &box_direction,
                            std::string box_position,
                            double m_dQdxRectangleLength,
                            double m_dQdxRectangleWidth,
                            std::string _pfp_producer)
{
  if (box_position != "start" && box_position != "end")
  {
    return;
  }

  std::vector<double> _gain;
  if (evt.isRealData())
    _gain = _data_gain;
  else
    _gain = _mc_gain;

  detinfo::DetectorProperties const *detprop =
      lar::providerFrom<detinfo::DetectorPropertiesService>();
  // art::ServiceHandle<geo::Geometry> geom;

  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);

  TVector3 pfp_dir;

  //For a shower
  if (pfparticle_handle->at(pfp_id).PdgCode() == 11)
  {
    try
    {
      art::FindOneP<recob::Shower> shower_per_pfpart(pfparticle_handle, evt, _pfp_producer);
      auto const &shower_obj = shower_per_pfpart.at(pfp_id);
      pfp_dir.SetX(shower_obj->Direction().X());
      pfp_dir.SetY(shower_obj->Direction().Y());
      pfp_dir.SetZ(shower_obj->Direction().Z());
    }
    catch (...)
    {
      return;
    }
  }
  // For a track
  else
  {
    try
    {
      art::FindOneP<recob::Track> track_per_pfpart(pfparticle_handle, evt, _pfp_producer);
      auto const &track_obj = track_per_pfpart.at(pfp_id);

      pfp_dir.SetX(track_obj->StartDirection().X());
      pfp_dir.SetY(track_obj->StartDirection().Y());
      pfp_dir.SetZ(track_obj->StartDirection().Z());
    }
    catch (...)
    {
      return;
    }
  }

  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt,
                                                     _pfp_producer);
  art::FindManyP<recob::Hit> hits_per_clusters(cluster_handle, evt,
                                               _pfp_producer);

  std::vector<art::Ptr<recob::Cluster>> clusters = clusters_per_pfpart.at(pfp_id);

  for (size_t icl = 0; icl < clusters.size(); icl++)
  {
    std::vector<art::Ptr<recob::Hit>> hits =
        hits_per_clusters.at(clusters[icl].key());

    double wireSpacing = 0.3;
    double tolerance = 0.001;
    std::vector<double> cluster_axis;
    std::vector<double> cluster_start;
    std::vector<double> cluster_end;

    double pitch = geoHelper.getPitch(pfp_dir, clusters[icl]->Plane().Plane);
    if (box_position == "end")
    {
      pitch *= -1;
    }

    double start_x = detprop->ConvertTicksToX(clusters[icl]->StartTick(), clusters[icl]->Plane());
    double end_x = detprop->ConvertTicksToX(clusters[icl]->EndTick(), clusters[icl]->Plane());
    // std::cout << "start " << start_x << ", end: " << end_x << std::endl;
    if (pitch >= 0)
    {
      std::reverse(hits.begin(), hits.end());
      cluster_axis = {cos(clusters[icl]->StartAngle()),
                      sin(clusters[icl]->StartAngle())};
      cluster_start = {clusters[icl]->StartWire() * wireSpacing - tolerance * cos(clusters[icl]->StartAngle()),
                       start_x - tolerance * sin(clusters[icl]->StartAngle())};
      cluster_end = {clusters[icl]->EndWire() * wireSpacing,
                     end_x};
    }
    else
    {
      cluster_axis = {-1. * cos(clusters[icl]->StartAngle()),
                      -1. * sin(clusters[icl]->StartAngle())};
      cluster_start = {clusters[icl]->EndWire() * wireSpacing + tolerance * cos(clusters[icl]->StartAngle()),
                       end_x + tolerance * sin(clusters[icl]->StartAngle())};
      cluster_end = {clusters[icl]->StartWire() * wireSpacing,
                     start_x};
    }

    double cluster_length = sqrt(pow(cluster_end[0] - cluster_start[0], 2) +
                                 pow(cluster_end[1] - cluster_start[1], 2));
    if (cluster_length <= 0)
      continue;

    // Build rectangle 4 x 1 cm around the cluster axis
    std::vector<std::vector<double>> points;
    geoHelper.buildRectangle(m_dQdxRectangleLength, m_dQdxRectangleWidth,
                             cluster_start, cluster_axis, points);

    //std::cout << "[dQdx] Point 1 " << points[0][0] << " " << points[0][1] << std::endl;
    //std::cout << "[dQdx] Point 2 " << points[1][0] << " " << points[1][1] << std::endl;
    //std::cout << "[dQdx] Point 3 " << points[2][0] << " " << points[2][1] << std::endl;
    //std::cout << "[dQdx] Point 4 " << points[3][0] << " " << points[3][1] << std::endl;

    std::vector<double> dqdxs;

    for (auto &hit : hits)
    {
      std::vector<double> hit_pos = {hit->WireID().Wire * wireSpacing,
                                     detprop->ConvertTicksToX(hit->PeakTime(), clusters[icl]->Plane())};
      bool is_within = geoHelper.isInside(hit_pos, points);

      if (is_within)
      {
        double q = integral2charge(evt, hit->Integral(), clusters[icl]->Plane().Plane);
        double dedx = dQdx2dEdx(q / fabs(pitch));
        dqdxs.push_back(dedx);
        if (clusters[icl]->Plane().Plane == 2)
        {
          dqdx_hits.push_back(dedx);
          dqdx_wires.push_back(hit->WireID().Wire);
          box_start[0] = cluster_start[0];
          box_start[1] = cluster_start[1];
          box_direction[0] = cluster_axis[0];
          box_direction[1] = cluster_axis[1];
        }
      }
    }

    // Get the median
    if (dqdxs.size() > 0)
    {
      std::nth_element(dqdxs.begin(), dqdxs.begin() + dqdxs.size() / 2, dqdxs.end());
      dqdx[clusters[icl]->Plane().Plane] = dqdxs[dqdxs.size() / 2];
    }
  }
}

void EnergyHelper::dEdxFromdQdx(std::vector<double> &dedx,
                                std::vector<double> &dqdx)
{
  for (size_t i = 0; i < dqdx.size(); i++)
  {
    if (dqdx[i] > 0)
      dedx[i] = dqdx[i] * (work_function) / recombination_factor;
  }
}

double EnergyHelper::integral2charge(const art::Event &evt, double const integral, int const plane) const
{
  if (evt.isRealData())
  {
    return (integral * _data_gain[plane]);
  }
  else
  {
    return (integral * _mc_gain[plane]);
  }
}

double EnergyHelper::dQdx2dEdx(double const dqdx) const
{
  return (dqdx * work_function / recombination_factor);
}

} // namespace lee

#endif