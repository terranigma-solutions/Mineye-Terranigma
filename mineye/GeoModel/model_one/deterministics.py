from dataclasses import dataclass

@dataclass(frozen=True)
class DetGravityResponseAligned:
    norm_params: dict
    def __call__(self, samples, geo_model, solutions):
        return align_forward_to_observed(-solutions.gravity, self.norm_params)

@dataclass(frozen=True)
class DetMeanGravityAligned:
    norm_params: dict
    def __call__(self, samples, geo_model, solutions):
        mu = align_forward_to_observed(-solutions.gravity, self.norm_params)
        return torch.mean(mu)

@dataclass(frozen=True)
class DetMaxGravityAligned:
    norm_params: dict
    def __call__(self, samples, geo_model, solutions):
        mu = align_forward_to_observed(-solutions.gravity, self.norm_params)
        # torch.max(..., dim=0) returns (values, indices) — pyro.deterministic needs a Tensor.
        values, idx = torch.max(mu, dim=0)
        # If you want both, register two deterministics (see dict below).
        return values

@dataclass(frozen=True)
class DetArgmaxGravityAligned:
    norm_params: dict
    def __call__(self, samples, geo_model, solutions):
        mu = align_forward_to_observed(-solutions.gravity, self.norm_params)
        _, idx = torch.max(mu, dim=0)
        return idx