# BLAKE2b — RFC 7693
# Obtained from: https://github.com/todoforai/BLAKE2.jl under MIT License
# Git hash: 096b74797207c1d9e789639d3ade8482ff5242ec

module Blake2b

    using StaticArrays: MVector, SVector

    const BLAKE2B_BLOCKBYTES = 128
    const BLAKE2B_OUTBYTES   = 64
    const BLAKE2B_KEYBYTES   = 64

    const BLAKE2B_IV = (
        0x6a09e667f3bcc908, 0xbb67ae8584caa73b,
        0x3c6ef372fe94f82b, 0xa54ff53a5f1d36f1,
        0x510e527fade682d1, 0x9b05688c2b3e6c1f,
        0x1f83d9abfb41bd6b, 0x5be0cd19137e2179,
    )

    const BLAKE2B_SIGMA = (
        ( 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15),
        (14, 10,  4,  8,  9, 15, 13,  6,  1, 12,  0,  2, 11,  7,  5,  3),
        (11,  8, 12,  0,  5,  2, 15, 13, 10, 14,  3,  6,  7,  1,  9,  4),
        ( 7,  9,  3,  1, 13, 12, 11, 14,  2,  6,  5, 10,  4,  0, 15,  8),
        ( 9,  0,  5,  7,  2,  4, 10, 15, 14,  1, 11, 12,  6,  8,  3, 13),
        ( 2, 12,  6, 10,  0, 11,  8,  3,  4, 13,  7,  5, 15, 14,  1,  9),
        (12,  5,  1, 15, 14, 13,  4, 10,  0,  7,  6,  3,  9,  2,  8, 11),
        (13, 11,  7, 14, 12,  1,  3,  9,  5,  0, 15,  4,  8,  6,  2, 10),
        ( 6, 15, 14,  9, 11,  3,  0,  8, 12,  2, 13,  7,  1,  4, 10,  5),
        (10,  2,  8,  4,  7,  6,  1,  5, 15, 11,  9, 14,  3, 12, 13,  0),
        ( 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15),
        (14, 10,  4,  8,  9, 15, 13,  6,  1, 12,  0,  2, 11,  7,  5,  3),
    )

    mutable struct Blake2bContext
        h::MVector{8, UInt64}
        t::UInt128              # byte counter
        f::UInt64               # finalization flag
        buf::MVector{128, UInt8}
        buflen::Int
        outlen::Int
    end

    """
        Blake2bContext(; output_length=64, key=UInt8[])

    Create an incremental BLAKE2b hashing context.
    """
    function Blake2bContext(; output_length::Int=BLAKE2B_OUTBYTES, key::AbstractVector{UInt8}=UInt8[])
        1 <= output_length <= BLAKE2B_OUTBYTES || throw(ArgumentError("output_length must be 1–$BLAKE2B_OUTBYTES"))
        length(key) <= BLAKE2B_KEYBYTES || throw(ArgumentError("key must be ≤ $BLAKE2B_KEYBYTES bytes"))

        h = MVector{8, UInt64}(BLAKE2B_IV...)
        h[1] ⊻= UInt64(output_length) | (UInt64(length(key)) << 8) | 0x0000000001010000

        buf = zero(MVector{128, UInt8})
        ctx = Blake2bContext(h, UInt128(0), UInt64(0), buf, 0, output_length)

        if !isempty(key)
            block = zero(MVector{128, UInt8})
            @inbounds for i in 1:length(key); block[i] = key[i]; end
            update!(ctx, block)
        end
        ctx
    end

    @inline rotr64(x::UInt64, n::Int) = (x >> n) | (x << (64 - n))

    @inline function _gb!(v, a, b, c, d, x::UInt64, y::UInt64)
        v[a] = v[a] + v[b] + x;  v[d] = rotr64(v[d] ⊻ v[a], 32)
        v[c] = v[c] + v[d];      v[b] = rotr64(v[b] ⊻ v[c], 24)
        v[a] = v[a] + v[b] + y;  v[d] = rotr64(v[d] ⊻ v[a], 16)
        v[c] = v[c] + v[d];      v[b] = rotr64(v[b] ⊻ v[c], 63)
        nothing
    end

    function _blake2b_compress!(ctx::Blake2bContext, block::AbstractVector{UInt8})
        m = MVector{16, UInt64}(undef)
        @inbounds for i in 1:16
            off = (i - 1) * 8
            m[i] = UInt64(block[off+1])       | (UInt64(block[off+2]) << 8)  |
                (UInt64(block[off+3]) << 16) | (UInt64(block[off+4]) << 24) |
                (UInt64(block[off+5]) << 32) | (UInt64(block[off+6]) << 40) |
                (UInt64(block[off+7]) << 48) | (UInt64(block[off+8]) << 56)
        end

        v = MVector{16, UInt64}(undef)
        @inbounds for i in 1:8; v[i] = ctx.h[i]; end
        v[9]  = BLAKE2B_IV[1]; v[10] = BLAKE2B_IV[2]
        v[11] = BLAKE2B_IV[3]; v[12] = BLAKE2B_IV[4]
        v[13] = BLAKE2B_IV[5]; v[14] = BLAKE2B_IV[6]
        v[15] = BLAKE2B_IV[7]; v[16] = BLAKE2B_IV[8]

        v[13] ⊻= UInt64(ctx.t & 0xFFFFFFFFFFFFFFFF)
        v[14] ⊻= UInt64((ctx.t >> 64) & 0xFFFFFFFFFFFFFFFF)
        v[15] ⊻= ctx.f

        @inbounds for round in 1:12
            s = BLAKE2B_SIGMA[round]
            _gb!(v, 1, 5,  9, 13, m[s[ 1]+1], m[s[ 2]+1])
            _gb!(v, 2, 6, 10, 14, m[s[ 3]+1], m[s[ 4]+1])
            _gb!(v, 3, 7, 11, 15, m[s[ 5]+1], m[s[ 6]+1])
            _gb!(v, 4, 8, 12, 16, m[s[ 7]+1], m[s[ 8]+1])
            _gb!(v, 1, 6, 11, 16, m[s[ 9]+1], m[s[10]+1])
            _gb!(v, 2, 7, 12, 13, m[s[11]+1], m[s[12]+1])
            _gb!(v, 3, 8,  9, 14, m[s[13]+1], m[s[14]+1])
            _gb!(v, 4, 5, 10, 15, m[s[15]+1], m[s[16]+1])
        end

        @inbounds for i in 1:8
            ctx.h[i] ⊻= v[i] ⊻ v[i+8]
        end
        nothing
    end

    """
        update!(ctx::Blake2bContext, data::AbstractVector{UInt8})

    Feed data into the BLAKE2b context.
    """
    function update!(ctx::Blake2bContext, data::AbstractVector{UInt8})
        off = 1
        remaining = length(data)
        while remaining > 0
            if ctx.buflen == BLAKE2B_BLOCKBYTES
                ctx.t += BLAKE2B_BLOCKBYTES
                _blake2b_compress!(ctx, ctx.buf)
                ctx.buflen = 0
            end
            n = min(remaining, BLAKE2B_BLOCKBYTES - ctx.buflen)
            @inbounds for i in 0:n-1; ctx.buf[ctx.buflen + 1 + i] = data[off + i]; end
            ctx.buflen += n
            off += n
            remaining -= n
        end
        ctx
    end

    # ─── Static-friendly core: writes into caller's buffer ───────────────────────

    """
        digest!(out::AbstractVector{UInt8}, ctx::Blake2bContext) -> Nothing

    Finalize and write the hash into `out`. The context should not be reused after this.
    `out` must have at least `ctx.outlen` bytes.
    """
    function digest!(out::AbstractVector{UInt8}, ctx::Blake2bContext)
        ctx.t += ctx.buflen
        ctx.f = 0xFFFFFFFFFFFFFFFF
        @inbounds for i in (ctx.buflen + 1):BLAKE2B_BLOCKBYTES
            ctx.buf[i] = 0x00
        end
        _blake2b_compress!(ctx, ctx.buf)

        @inbounds for i in 1:((ctx.outlen + 7) ÷ 8)
            off = (i - 1) * 8
            w = ctx.h[i]
            for j in 1:min(8, ctx.outlen - off)
                out[off+j] = UInt8((w >> (8*(j-1))) & 0xFF)
            end
        end
        nothing
    end

    # ─── Convenience wrappers (allocating) ───────────────────────────────────────

    """
        digest!(ctx::Blake2bContext) -> Vector{UInt8}

    Finalize and return the hash. The context should not be reused after this.
    """
    function digest!(ctx::Blake2bContext)
        out = Vector{UInt8}(undef, ctx.outlen)
        digest!(out, ctx)
        out
    end

    """
        compute2b(data; output_length=64, key=UInt8[]) -> Vector{UInt8}

    Compute BLAKE2b hash of `data`.

    # Arguments
    - `output_length::Int=64`: Output hash length in bytes (1–64)
    - `key::AbstractVector{UInt8}=UInt8[]`: Optional key for keyed hashing (MAC)
    """
    function compute2b(data::AbstractVector{UInt8};
                    output_length::Int=BLAKE2B_OUTBYTES,
                    key::AbstractVector{UInt8}=UInt8[])
        ctx = Blake2bContext(; output_length, key)
        update!(ctx, data)
        digest!(ctx)
    end
    export compute2b

end # end module