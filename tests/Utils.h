#pragma once

#include <Hull/Hull.h>
#include <memory>
#include <string>

namespace hull {
void toObj(const std::vector<Coordinate> &vertices,
           const std::vector<Facet> &facets, const std::string &fileName);

bool check_normals(const hull::Hull &subject);

class TrivialObserver : public Observer {
public:
  TrivialObserver() = default;

  std::unique_ptr<Notification> last_notification;
  void hullChanges(const Notification &notification) override {
    last_notification = std::make_unique<Notification>(notification);
  };
};
} // namespace hull

std::string generate_obj_log_name(const std::string &file_name);
