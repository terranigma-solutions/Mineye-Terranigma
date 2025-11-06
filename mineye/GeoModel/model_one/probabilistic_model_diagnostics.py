import io

import pyro
from PIL import Image
from matplotlib import pyplot
from pyro.infer.inspect import get_dependencies


def debug_gradient_flow(likelihood_fn, solutions, obs_data):
    """
    Systematically check where gradient flow breaks.
    
    Call this BEFORE running NUTS to diagnose issues.
    """
    import torch

    print("\n" + "=" * 80)
    print("GRADIENT FLOW DIAGNOSTIC")
    print("=" * 80)

    # Get the distribution
    lik_dist = likelihood_fn(solutions)

    # Compute log probability
    try:
        log_prob = lik_dist.log_prob(obs_data)
        print(f"✓ Log probability computed: {log_prob.item():.4f}")
    except Exception as e:
        print(f"✗ Failed to compute log_prob: {e}")
        return

    # Check if we can compute gradients
    print(f"\nLog prob requires_grad: {log_prob.requires_grad}")
    print(f"Log prob grad_fn: {log_prob.grad_fn}")

    if not log_prob.requires_grad:
        print("\n⚠️  WARNING: log_prob does not require gradients!")
        print("This means NUTS cannot compute gradients for sampling.")

        # Trace back through the computation graph
        print("\nTracing back through computation:")
        print(f"  Distribution type: {type(lik_dist).__name__}")

        # Check distribution parameters
        if hasattr(lik_dist, 'loc'):
            print(f"  loc requires_grad: {lik_dist.loc.requires_grad}")
            print(f"  loc grad_fn: {lik_dist.loc.grad_fn}")

        if hasattr(lik_dist, 'scale_tril'):
            print(f"  scale_tril requires_grad: {lik_dist.scale_tril.requires_grad}")
            print(f"  scale_tril grad_fn: {lik_dist.scale_tril.grad_fn}")

        if hasattr(lik_dist, 'covariance_matrix'):
            print(f"  covariance_matrix requires_grad: {lik_dist.covariance_matrix.requires_grad}")
            print(f"  covariance_matrix grad_fn: {lik_dist.covariance_matrix.grad_fn}")

        if hasattr(lik_dist, 'df'):
            print(f"  df requires_grad: {lik_dist.df.requires_grad}")
            print(f"  df grad_fn: {lik_dist.df.grad_fn}")
    else:
        print("\n✓ Gradient flow is intact!")

        # Try computing gradients
        try:
            # Create dummy parameters to test backward pass
            test_grad = torch.autograd.grad(
                outputs=log_prob,
                inputs=[p for p in lik_dist.parameters() if p.requires_grad],
                allow_unused=True
            )
            print(f"✓ Backward pass successful, got {len(test_grad)} gradients")
        except Exception as e:
            print(f"✗ Backward pass failed: {e}")

    print("=" * 80 + "\n")


def trace_pyro_model(prob_model, geo_model, obs_data, print_full=True):
    """
    Use Pyro's trace functionality to see all sample sites and deterministics.

    This shows you EXACTLY what Pyro sees.
    """
    from pyro import poutine
    import torch

    print("\n" + "=" * 80)
    print("PYRO MODEL TRACE")
    print("=" * 80 + "\n")

    # Trace the model execution
    trace = poutine.trace(prob_model).get_trace(geo_model, obs_data)

    print(f"{'Site Name':<30} | {'Type':<15} | {'Shape':<15} | {'Grad?':<10} | {'Grad Fn'}")
    print("-" * 100)

    for name, node in trace.nodes.items():
        if node['type'] == 'sample':
            value = node['value']

            if isinstance(value, torch.Tensor):
                shape_str = str(tuple(value.shape))
                grad_str = "✓" if value.requires_grad else "✗"
                grad_fn_str = str(value.grad_fn)[:30] if value.grad_fn else "None"

                print(f"{name:<30} | {node['type']:<15} | {shape_str:<15} | {grad_str:<10} | {grad_fn_str}")

                # Flag problematic ones
                if not value.requires_grad and name not in ['obs', 'Gravity Measurement']:
                    print(f"  ⚠️  WARNING: Sample site '{name}' has no gradient flow!")
            else:
                print(f"{name:<30} | {node['type']:<15} | {str(type(value)):<15} | {'N/A':<10} | N/A")

    if print_full:
        print("\n" + "=" * 80)
        print("FULL TRACE DETAILS")
        print("=" * 80 + "\n")
        print(trace.format_shapes())

    print("\n" + "=" * 80 + "\n")

    return trace


def _plot_probability_model_graph(model, geo_model, y_obs_list=None):
    # ! This is not working well, the geo_model dependency is not properly picked up
    # ! This does not work at all when keops is on
    
    from pyro.infer.inspect import get_dependencies
    import pyro
    dependencies = get_dependencies(model, model_args=(geo_model, y_obs_list[:1]))
    dependencies

    # %%
    graph = pyro.render_model(
        model=model,
        model_args=(geo_model, y_obs_list,),
        render_params=True,
        render_distributions=True,
        render_deterministic=True
    )

    graph.attr(dpi='300')
    # Convert the graph to a PNG image format
    s = graph.pipe(format='png')

    # Open the image with PIL
    from PIL import Image
    image = Image.open(io.BytesIO(s))

    # Plot the image with matplotlib
    plt.figure(figsize=(10, 4))
    plt.imshow(image)
    plt.axis('off')  # Turn off axis
    plt.show()