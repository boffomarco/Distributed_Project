function [ p, t ] = cvt_circle_uniform ( n, sample_num, it_num, ...
  delaunay_display )

%*****************************************************************************80
%
%% CVT_CIRCLE_UNIFORM computes a CVT in a circle with uniform density.
%
%  Discussion:
%
%    This simple example carries out an iterative CVT calculation in a circle
%    with a uniform sampling density. 
%
%    This example has been modified to respond to MATLAB's undesired
%    changes in the Delaunay and Voronoi calculation procedures.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    20 July 2016
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Qiang Du, Vance Faber, Max Gunzburger,
%    Centroidal Voronoi Tessellations: Applications and Algorithms,
%    SIAM Review, 
%    Volume 41, 1999, pages 637-676.
%
%  Parameters:
%
%    Input, integer N, the number of generators.
%
%    Input, integer SAMPLE_NUM, the number of sample points.
%
%    Input, integer IT_NUM, the number of iterations to take.
%
%    Input, logical DELAUNAY_DISPLAY, is true if the Delaunay
%    triangulation is to be displayed.
%
%    Output, real P(2,N), the location of the generators.
%
%    Output, integer T(NT,3), information defining the Delaunay
%    triangulation of the generators.  NT is the number of triangles,
%    which varies depending on the arrangement of the generators.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'CVT_CIRCLE_UNIFORM:\n' );
  fprintf ( 1, '  A simple demonstration of a CVT computation\n' );
  fprintf ( 1, '  (Centroidal Voronoi Tessellation)\n' );
  fprintf ( 1, '  in a circle, with a uniform density.\n' );
  fprintf ( 1, '\n' );

  if ( nargin < 1 )
    n = 100;
    fprintf ( 1, '  The default number of generators N = %d\n', n );
  else
    fprintf ( 1, '  The user number of generators N = %d\n',  n );
  end

  if ( nargin < 2 )
    sample_num = 10000 * n;
    fprintf ( 1, '  The default number of sample points SAMPLE_NUM = %d\n', ...
      sample_num );
  else
    fprintf ( 1, '  The user number of sample points SAMPLE_NUM = %d\n', ...
      sample_num );
  end

  if ( nargin < 3 )
    it_num = 20;
    fprintf ( 1, '  The default number of iterations IT_NUM = %d.\n', ...
      it_num );
  else
    fprintf ( 1, '  The user number of iterations = %d\n', ...
      it_num );
  end

  if ( nargin < 4 )
    delaunay_display = true;
    fprintf ( 1, '  The default value of DELAUNAY_DISPLAY = %d.\n', ...
      delaunay_display );
  else
    fprintf ( 1, '  The user value of DELAUNAY_DISPLAY = %d\n', ...
      delaunay_display );
  end
%
%  This switch is set to true if the ACCUMARRAY command is available.
%  That speeds up the calculation a lot.  If you don't have the ACCUMARRAY
%  command, just set this to false.
%  
  accumarray_available = true;

  if ( accumarray_available )
    fprintf ( 1, '\n' );
    fprintf ( 1, '  The ACCUMARRAY command will be called.\n' );
  else
    fprintf ( 1, '\n' );
    fprintf ( 1, '  The ACCUMARRAY command will NOT be called.\n' );
  end
%
%  Clear the figure screen, if already open.
%
  clf
%
%  Initialize the generators.
%
  p = circle_uniform ( n );

  it = 0;
  
  while ( it <= it_num )
%
%  Compute the Delaunay triangle information T for the current nodes.
%
    if ( true )
      t = delaunay ( p(:,1), p(:,2) );
    else
      t = DelaunayTri ( p(:,1), p(:,2) );
    end
%
%  Display the Delaunay triangulation, if requested.
%
    if ( delaunay_display )
      figure ( 1 )
      triplot ( t, p(:,1), p(:,2) );
      title_string = sprintf ( 'Delaunay, step %d', it );
      title ( title_string );
      axis equal
      view ( 2 )
    end
%
%  Display the CVT generators, and the associated Voronoi diagram.
%
    figure ( 2 )
    clf
    hold on
    theta = 2.0 * pi * ( 0 : 50 ) / 50;
    x = cos ( theta );
    y = sin ( theta );

    if ( true )
      voronoi ( p(:,1), p(:,2) );
%     voronoi ( p(:,1), p(:,2), t );
    else
      plot ( p(:,1), p(:,2), 'ko', 'Markersize', 10 );
    end

    line ( x, y );
    title_string = sprintf ( 'Uniform Density on a Circle, CVT step %d', it );
    title ( title_string );
    axis equal
    drawnow
    hold off
%
%  Generate sample points.  
%  
    ps = circle_uniform ( sample_num );
%
%  For each sample point, find K, the index of the nearest generator.
%  We do this efficiently by using the Delaunay information with
%  Matlab's DSEARCHN command, rather than a brute force nearest neighbor
%  computation.
%
    if ( true )
      k = dsearchn ( p, t, ps );
    else
      k = nearestNeighbor ( t, ps(:,1), ps(:,2) );
    end
%
%  The centroid of the Voronoi region associated with each generator
%  is approximated by the average of the sample points it was closest to.
%
    if ( accumarray_available )

      count(1:n,1) = accumarray ( k, ones(sample_num,1) );
      centroid(1:n,1) = accumarray ( k, ps(:,1) );
      centroid(1:n,2) = accumarray ( k, ps(:,2) );

    else

      count(1:n,1) = 0;
      centroid(1:n,1) = 0.0;
      centroid(1:n,2) = 0.0;

      for i = 1 : sample_num
        j = k(i);
        count(j,1) = count(j,1) + 1;
        centroid(j,1) = centroid(j,1) + ps(i,1);
        centroid(j,2) = centroid(j,2) + ps(i,2);
      end

    end
%
%  Replace the generators by the centroids.
%
    p(1:n,1) = centroid(1:n,1) ./ count(1:n,1);
    p(1:n,2) = centroid(1:n,2) ./ count(1:n,1);

    it = it + 1;
    
  end
%
%  Save the graphics.
%
  figure ( 1 );
  filename = 'delaunay.png';
  print ( '-dpng', filename );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );

  figure ( 1 );
  filename = 'voronoi.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'CVT_CIRCLE_UNIFORM:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );

  return
end
function p = circle_uniform ( n )

%*****************************************************************************80
%
%% CIRCLE_UNIFORM uniformly samples points within a circle.
%
%  Discussion:
%
%    This routine returns N points sampled uniformly at random
%    from within the unit circle.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    20 July 2016
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the number of points to generate.
%
%    Output, real P(N,2), the sample points.
%

%
%  Generate N normally distributed points in R^2.
%
  p = randn ( n, 2 );
%
%  Normalize them (now they lie on the circumference of the unit circle).
%
  p_norm(1:n,1) = sqrt ( dot ( p, p, 2 ) );
%
%  Choose a value R uniformly in [0,1].
%
  r(1:n,1) = rand ( n, 1 );
%
%  The radius is the square root of R.
%
  for i = 1 : 2
    p(1:n,i) = sqrt ( r(1:n) ) .* p(1:n,i) ./ p_norm(1:n);
  end

  return
end
function timestamp ( )

%*****************************************************************************80
%
%% TIMESTAMP prints the current YMDHMS date as a timestamp.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    14 February 2003
%
%  Author:
%
%    John Burkardt
%
  t = now;
  c = datevec ( t );
  s = datestr ( c, 0 );
  fprintf ( 1, '%s\n', s );

  return
end
