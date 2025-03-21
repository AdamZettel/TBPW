Subroutine fft(grid, nn, ndim, isign)

  use sysparams, only : double, twopi, zero, half, one, two, imaginary

  implicit none

  integer, intent(in)            :: isign, ndim, nn(ndim)
  complex(double), intent(inout) :: grid(:)

  integer i1, i2, i2rev, i3, i3rev, ibit, idim, ifp1, ifp2, ip1, ip2, &
        & ip3, k1, k2, n1, nprev, nrem, ntot
  complex(double) temp, w, wp, wtemp
  real(double) theta

  ntot = 1
  do idim =1, ndim
    ntot = ntot * nn(idim)
  end do

  nprev = 1

  do idim = 1, ndim

    n1 = nn(idim)
    nrem = ntot/(n1*nprev)
    ip1 = nprev
    ip2 = ip1 * n1
    ip3 = ip2 * nrem
    i2rev = 1

    do i2 = 1, ip2, ip1 
      if(i2 .lt. i2rev) then
        do i1 = i2, i2+ip1-1
          do i3 = i1, ip3, ip2

            i3rev = i2rev + i3 - i2
            temp = grid(i3)
            grid(i3) = grid(i3rev)
            grid(i3rev) = temp

          end do
        end do
      end if

      ibit = ip2/2

      do while((ibit .ge. ip1) .and. (i2rev .gt. ibit))

        i2rev = i2rev - ibit
        ibit = ibit / 2

      end do

      i2rev = i2rev + ibit

    end do

    ifp1 = ip1

    do while (ifp1 .lt. ip2)

      ifp2 = 2 * ifp1
      theta = isign * twopi /(ifp2 / ip1)
      wp = -two * sin(half * theta) ** 2 + imaginary * sin(theta)
      w = one
 
      do i3 = 1, ifp1, ip1
        do i1 = i3, i3 + ip1 - 1, 1
          do i2 = i1, ip3, ifp2

            k1 = i2
            k2 = k1 + ifp1
			temp = w * grid(k2)
            grid(k2) = grid(k1) - temp
            grid(k1) = grid(k1) + temp

          end do
        end do

        w = w + w * wp

      end do

      ifp1 = ifp2

    end do

    nprev = n1 * nprev

  end do

  return

end Subroutine fft
