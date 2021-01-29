open System
open System.IO
open Aardvark.Base

module Histogram =
    /// get a histogram for an image
    let compute (img : PixImage<uint16>) =
        let histogram = Array.zeroCreate 65536
        let channel = img.GetChannel(0L)
        channel.ForeachIndex(fun i -> histogram.[int channel.[i]] <- histogram.[int channel.[i]] + 1L) |> ignore
        histogram

    /// get a otsu-threshold based on a histogram (as a floating point value in [0..1])
    let getOtsuThreshold (histo : int64[]) =
        let mutable lsum = 0.0
        let mutable lsumSq = 0.0
        let mutable lw = 0.0
                
        let mutable rsum = 0.0
        let mutable rsumSq = 0.0
        let mutable rw = 0.0

        for i in 0 .. histo.Length - 1 do
            let v = float histo.[i]
            rsum <- rsum + float i * v
            rsumSq <- rsumSq + (float i * float i) * v
            rw <- rw + v
                    
                    
        let mutable bestWeight = System.Double.PositiveInfinity
        let mutable bestThreshold = -1
        for i in 0 .. histo.Length - 1 do
            let lv = 
                if lw <= 0.0 then 
                    0.0
                else
                    let avg = (lsum / lw)
                    let avgSq = (lsumSq / lw)
                    avgSq - avg*avg
                        
            let rv = 
                if rw <= 0.0 then
                    0.0
                else
                    let avg = (rsum / rw)
                    let avgSq = (rsumSq / rw)
                    avgSq - avg*avg

            let w = (lw * lv + rw * rv) / (lw + rw)
                    
            if w < bestWeight then
                bestWeight <- w
                bestThreshold <- i
        
            let v = float histo.[i]
            let a = float i * v
            let b = (float i * float i) * v
            lsum <- lsum + a
            lsumSq <- lsumSq + b
            lw <- lw + v
            rsum <- rsum - a
            rsumSq <- rsumSq - b
            rw <- rw - v
                    

        float bestThreshold / float histo.Length
        
module P202 =   
    // random generator
    let rand = RandomSystem()

    /// axis-aligned, quadratic, fixed-size, stochastic P202
    let growing (quadSize : int) (iterations : int) (binary : Matrix<byte>) =
        
        let inline sx (b : Box2i) = b.Max.X - b.Min.X + 1
        let inline sy (b : Box2i) = b.Max.Y - b.Min.Y + 1


        let rec grow (size : int) (pixel : V2i) (mat : NativeMatrix<byte>) =
            let mutable box = Box2i(pixel, pixel)
            let mutable inside = 1
            let mutable failed = false

            // repeat until box is large enough
            while not failed && (sx box < size || sy box < size) do

                let mutable alternatives = []

                // test -x
                if box.Min.X > 0 && sx box < size then
                    // can go to -x
                    let l = box.Min.X - 1
                    
                    // count "1" pixels in border
                    let mutable gain = 0
                    mat.[l, box.Min.Y .. box.Max.Y] |> NativeVector.iter (fun _ a -> if a > 127uy then gain <- gain + 1)

                    // add alternative
                    alternatives <- (gain, Box2i(box.Min - V2i.IO, box.Max)) :: alternatives

                
                // test +x
                if box.Max.X < int (mat.SX - 1L) && sx box < size then
                    // can go to +x
                    let r = box.Max.X + 1

                    // count "1" pixels in border
                    let mutable gain = 0
                    mat.[r, box.Min.Y .. box.Max.Y] |> NativeVector.iter (fun _ a -> if a > 127uy then gain <- gain + 1)
                
                    // add alternative
                    alternatives <- (gain, Box2i(box.Min, box.Max + V2i.IO)) :: alternatives

                
                // test -y
                if box.Min.Y > 0 && sy box < size then
                    // can go to -y
                    let b = box.Min.Y - 1
                    
                    // count "1" pixels in border
                    let mutable gain = 0
                    mat.[box.Min.X .. box.Max.X, b] |> NativeVector.iter (fun _ a -> if a > 127uy then gain <- gain + 1)
                
                    // add alternative
                    alternatives <- (gain, Box2i(box.Min - V2i.OI, box.Max)) :: alternatives
                   
                // test +y 
                if box.Max.Y < int (mat.SY - 1L) && sy box < size then
                    // can go to +y
                    let t = box.Max.Y + 1
                    
                    // count "1" pixels in border
                    let mutable gain = 0
                    mat.[box.Min.X .. box.Max.X, t] |> NativeVector.iter (fun _ a -> if a > 127uy then gain <- gain + 1)
                
                    // add alternative
                    alternatives <- (gain, Box2i(box.Min, box.Max + V2i.OI)) :: alternatives
                  
            

                if List.isEmpty alternatives then
                    // no alternatives => image too small
                    failed <- true
                else
                    let mutable bestGain = 0
                    for (gain, _) in alternatives do bestGain <- max bestGain gain

                    // find the (possibly multiple) alternatives with optimal gain
                    let optimal = alternatives |> List.filter (fun (gain,_) -> gain = bestGain)


                    match optimal with
                    | [(bestGain, newBestBox)] ->
                        // unique best alternative => take
                        // maintain inside and update box
                        inside <- inside + bestGain
                        box <- newBestBox

                    | optimal ->
                        // multiple best alternatives
                        // => randomly choose one
                        let _, randomBox = optimal.[rand.UniformInt optimal.Length]
                        box <- randomBox


            if not failed then
                Some (inside, box)
            else
                None

        let rand = RandomSystem()
        let box, count = 
            NativeMatrix.using binary (fun pBinary ->
                let mutable bestBox = Box2i.Invalid
                let mutable bestCnt = 0
                let s = V2i pBinary.Size

                let mutable iters = 0
                //Log.startTimed "growing"
                while iters < iterations do
                    let px = rand.UniformV2i s
                    if pBinary.[px] > 127uy then
                        match grow quadSize px pBinary with
                        | Some (cnt, box) -> 
                            if cnt > bestCnt then
                                //Log.line "%.2f%%" (100.0 * float cnt / float (sx box * sy box))
                                bestCnt <- cnt
                                bestBox <- box
                        | None ->
                            ()
                        iters <- iters + 1
                        Report.Progress(float iters / float iterations)
                //Log.stop()

                bestBox, bestCnt
            )

        let porosity = float count / float (sx box * sy box)

        box, porosity

    /// axis-aligned, quadratic, fixed-size brute-force P202
    let bruteForce (quadSize : int) (binary : Matrix<byte>) =
    
        let inline sx (b : Box2i) = b.Max.X - b.Min.X + 1
        let inline sy (b : Box2i) = b.Max.Y - b.Min.Y + 1
        let half = V2i(quadSize, quadSize) / 2
        let upperHalf = quadSize - V2i.II - half

        let l = half
        let u = V2i binary.Size - V2i.II - upperHalf

        NativeMatrix.using binary (fun pBinary ->
            let mutable bestCount = 0
            let mutable bestBox = Box2i.Invalid

            let total = (1 + u.X - l.X) * (1 + u.Y - l.Y)
            let mutable i = 0
            //Log.startTimed "brute force"
            for x in l.X .. u.X do
                for y in l.Y .. u.Y do
                    let c = V2i(x,y)
                    let b = Box2i(c - half, c + upperHalf)
                    let quad = pBinary.[b.Min.X .. b.Max.X, b.Min.Y .. b.Max.Y]
                    let mutable cnt = 0
                    quad |> NativeMatrix.iter (fun _ v -> if v > 127uy then cnt <- cnt + 1)
                    if cnt > bestCount then
                        bestCount <- cnt
                        bestBox <- b

                    i <- i + 1
                    Report.Progress(float i / float total)
            //Log.stop()

            let porosity = float bestCount / float (sx bestBox * sy bestBox)
            bestBox, porosity

        )



type Args =
    {
        inputFile : string
        outputFile : string
        iterations : int
        quadSize : int
        skipBruteForce : bool
    }

// argument parser
module Args =
    
    let inline (|ExistingPath|_|) (v : string) =
        if File.Exists v then Some v
        else None
        
    let inline (|FilePath|_|) (v : string) =
        if not (v.StartsWith "-") then Some v
        else None

    let inline (|Int|_|) (v : string) =
        match System.Int32.TryParse v with
        | (true, v) -> Some v
        | _ -> None

    let rec parse (args : list<string>) (a : Args) =
        match args with
        | [] -> Success a
        | "-i" :: Int iterations :: rest -> parse rest { a with iterations = iterations }
        | "-s" :: Int size :: rest -> parse rest { a with quadSize = size }
        | "-o" :: FilePath o :: rest -> parse rest { a with outputFile = o }
        | "-sb" :: rest -> parse rest { a with skipBruteForce = true }
        | FilePath input :: rest -> parse rest { a with inputFile = input }
        | a :: rest -> Error (sprintf "unknown argument %A" a)

[<EntryPoint;STAThread>]
let main argv = 
    Aardvark.Init()
    let defaultArgs = { inputFile = null; outputFile = null; iterations = 1000; quadSize = 355; skipBruteForce = false }
    match Args.parse (Array.toList argv) defaultArgs with
    | Error err ->
        Log.error "%s" err
        -1
    | Success args when isNull args.inputFile || not (File.Exists args.inputFile) ->
        Log.error "no input file given"
        -1
    | Success args when isNull args.outputFile ->
        Log.error "no output file given (specify via -o <path>)"
        -1
    | Success args ->

        Log.start "info"
        Log.line "input:      \"%s\"" args.inputFile
        Log.line "output:     \"%s\"" args.outputFile
        Log.line "iterations: %d" args.iterations
        Log.line "quadSize:   %d" args.quadSize
        Log.stop()

        Log.startTimed "reading image"
        // read the image, compute histogram and determine otsu-threshold
        let slice = 
            let img = PixImage.Create(args.inputFile)
            if img.Format = Col.Format.Gray then img.ToPixImage<uint16>(Col.Format.Gray)
            else 
                let dst = PixImage<uint16>(Col.Format.Gray, img.Size)
                dst.GetChannel(0L).Set(img.ToPixImage<uint16>().GetChannel(0L)) |> ignore
                dst

        let histogram = Histogram.compute slice
        let threshold = Histogram.getOtsuThreshold histogram * 65535.0 |> round |> uint16
        Log.line "size: %dx%d" slice.Size.X slice.Size.Y
        Log.line "otsu: %d" threshold
        Log.stop()

        // create a binary image using the otsu-threshold (stored in 8bit)
        let binary = slice.GetChannel(0L).Map(fun v -> if v >= threshold then 255uy else 0uy)

        // growing approach
        Log.startTimed "growing"
        let growingBox, growingPorosity = P202.growing args.quadSize args.iterations binary
        Log.line "%A" growingBox
        Log.line "%.2f%%" (100.0 * growingPorosity)
        Log.stop()

        let bruteForceBox = 
            if not args.skipBruteForce then
                // brute force approach (for reference)
                Log.startTimed "brute force"
                // known values for test image (takes about 4min to compute)
                //let bruteForceBox = Box2i(V2i(60, 632), V2i(414, 986))
                //let bruteForcePorosity = 0.1825
                let bruteForceBox, bruteForcePorosity = P202.bruteForce args.quadSize binary
                Log.line "%A" bruteForceBox
                Log.line "%.2f%%" (100.0 * bruteForcePorosity)
                Log.stop()
                Some bruteForceBox
            else
                None

        Log.start "writing output"
        // draw and save the result
        let res = PixImage<byte>(Col.Format.RGBA, V2i binary.Size)
        let m = res.GetMatrix<C4b>()
        m.SetMap(binary, fun v -> C4b(v,v,v,255uy)) |> ignore
        m.SetRectangle(growingBox.Min, growingBox.Max, C4b.Green)
        match bruteForceBox with
        | Some bruteForceBox -> m.SetRectangle(bruteForceBox.Min, bruteForceBox.Max, C4b.Red)
        | _ -> ()
        res.SaveAsImage args.outputFile
        Log.stop()
        0
