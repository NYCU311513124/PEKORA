
%{
sora = 1
輸入阻抗--------------------------負載阻抗
        |            |
        |            |
        |            |

sora = 2
輸入阻抗-----------------------------負載阻抗
        |            |            |
        |            |            |
        |            |            |
%}

tic

% 設定基本規格
meow = 1; % 計算的結構 : 1.CRLH phase shifter 2.阻抗匹配結構
if(meow==1)
    kanata = 9:0.01:11; % 總頻率範圍 (GHz)
    polka = zeros(1,length(kanata));
    polka(:) = 75; % 負載阻抗 (ohm)
    otonose = [9.8 10 10.2;35 0 -35]; % 最低頻率、中心頻率和最高頻率點的相位變化量 (degree)
    tokino = 2; % 是否要求中心頻率符合指定相位 : 1.是 2.否
    hoshimachi = 2; % 排序計算結果的優先順序 : 1.反射係數 2.相位變化 3.反射係數與相位差總和
else
    S1='S Parameter.csv'; % 匯入HFSS的S11
    kanata = readmatrix(S1,'Range','A2:A202'); % 讀取頻率範圍 (GHz)
    botan = readmatrix(S1,'Range','B2:B202'); % 讀取S11大小 (magnitude)
    lamy = readmatrix(S1,'Range','C2:C202'); % 讀取S11角度 (degree)
    nene = botan.* exp(1i*lamy/180*pi); % 轉換S11為複數表示式
    polka = 75 * (1+nene) ./ (1-nene); % 將S11轉換為天線輸入阻抗(負載阻抗)
    otonose = 0; % 不考慮相位變化
end
towa = 75; % 輸入阻抗 (ohm)
luna = [9.8 10.2]; % 計算頻率範圍 (GHz)
marine = -10; % 設計目標的反射係數 S11 (dB)
hime = 10; % 中心頻率 (GHz)

% 設定結構參數的計算數值範圍
bae = [70 130]; % 傳輸線的特性阻抗 (ohm)
fauna = [90 180]; % 傳輸線的電氣長度 (degree)
kronii = [70 130]; % stub的特性阻抗 (ohm)
mumei = [90 180]; % stub的電氣長度 (degree)

% 控制是否分析計算傳輸線或stub的特性阻抗與電氣長度
irys = zeros(4,2);
irys(:,1) = [1 1 1 1]; % 是否固定數值:0為單一數值，1為設定之範圍
irys(:,2) = [126 100 76 100]; % 固定的單一數值

% 設定計算控制項
sora = 1; % 電路結構 : 1.T型 2.π型
azki = 4; % 電路階數 (1.2.3.4...)
roboco = 2; % 選擇權重大小 : 1.小 2.大
miko = 2; % 計算結果是否以滿足阻抗匹配的頻率比例優先排序 : 1.否 2.是
suisei = 100; % 輸出結果的粒子數量

shiori = 100; % 重新計算的次數
iofi = 10; % particle filter迭代的次數
chino = 100; % 初始粒子的數量
bijou = 3; % 迭代時權重最大的粒子保存數量

gigi = 1; % 初始粒子分布情形 : 1.均勻分布 2.常態分布
raora = 2; % 是否依照權重大小重新取樣 1.否 2.是
cecilia = 2; % 粒子的運動狀態 : 1.隨機改變粒子位置 2.依照目前位置添加雜訊
liz = 4; % 重新取樣的方法 : 1.multinomial resampling 2.residual resampling 3.systematic resampling 4.stratified resampling
nodoka = 0.3; % 常態分布的標準差
su = 5; % 粒子數值添加雜訊的最大量
zeta = 2; % 是否重新取樣的條件 (nerissa/zeta)

myu = 2; % 是否執行KLD取樣
vitte = 0.1; % 與真實機率分布之間的誤差
rara = 1; % 方格的大小
nina = 2; % 是否對標準常態分布做修正 1.否 2.是
mugenshoujo = 1; % 標準常態分布的控制項

% 初始化計算參數
[watame,hina,todoroki,ui,oozora,ichijou] = initialize1(kanata,polka,luna,otonose,marine,hime,meow);
[okayu,haachama,ina,laplas,koyori,darknesss,hakui,kikirara] = initialize2(sora,azki,suisei,watame,irys,bae,fauna,kronii,mumei);

for k = 1:shiori % 每次重新計算particle filter
    nerissa = chino; % 粒子的數量
    if(meow==1) % CRLH phase shifter
        moona = zeros(4,nerissa); % 初始化每個粒子的參數
        if(gigi == 1) % 隨機產生均勻分布的粒子
            for i = 1:nerissa % 生成粒子
                for j = 1:4 % 隨機生成粒子數值
                    moona(j,i) = haachama(j,1) + (haachama(j,2)-haachama(j,1)) * rand;
                end
            end
        else % 隨機產生常態分布的粒子
            for i = 1:nerissa % 生成粒子
                for j = 1:4 % 隨機生成粒子數值
                    moona(j,i) = haachama(j,1) + (haachama(j,2)-haachama(j,1)) * (nodoka*randn+0.5);
                    if(moona(j,i) < haachama(j,1)) % 避免數值超出設定的範圍
                        moona(j,i) = haachama(j,1);
                    elseif(moona(j,i) > haachama(j,2))
                        moona(j,i) = haachama(j,2);
                    end
                end
            end
        end
    else % 阻抗匹配結構
        if(sora == 1) % T型電路結構
            moona =  zeros(4,okayu,nerissa); % 初始化每個粒子的參數
        else % π型電路結構
            moona =  zeros(4,azki-okayu,nerissa);
        end
        if(gigi == 1) % 隨機產生均勻分布的粒子
            for i = 1:nerissa % 生成粒子
                for j = 1:4 % 隨機生成粒子數值
                    if(j < 3) % 傳輸線參數
                        noel = okayu;
                    else % stub參數
                        noel = azki-okayu;
                    end
                    for m = 1:noel
                        moona(j,m,i) = haachama(j,1) + (haachama(j,2)-haachama(j,1)) * rand;
                    end
                end
            end
        else % 隨機產生常態分布的粒子
            for i = 1:nerissa % 生成粒子
                for j = 1:4 % 隨機生成粒子數值
                    if(j < 3) % 傳輸線參數
                        noel = okayu;
                    else % stub參數
                        noel = azki-okayu;
                    end
                    for m = 1:noel
                        moona(j,m,i) = haachama(j,1) + (haachama(j,2)-haachama(j,1)) * (nodoka*randn+0.5);
                        if(moona(j,m,i) < haachama(j,1)) % 避免數值超出設定的範圍
                            moona(j,m,i) = haachama(j,1);
                        elseif(moona(j,m,i) > haachama(j,2))
                            moona(j,m,i) = haachama(j,2);
                        end
                    end
                end
            end
        end
    end

    for h = 1:iofi % 每次particle filter迭代計算
        risu = zeros(1,nerissa); % 初始化每個粒子的反射係數 S11
        reine = zeros(1,nerissa); % 初始化每個粒子達到阻抗匹配的情形
        if(meow==1) % CRLH phase shifter
            isaki = zeros(1,nerissa); % 初始化每個粒子的相位差
            mizumiya = zeros(1,nerissa); % 初始化每個粒子的反射係數與相位差總和
            FLOWGLOW = zeros(1,nerissa); % 初始化每個粒子的相位變化斜率
        end
        
        for g = 1:nerissa % 對每個粒子做計算
            % 執行PEKO演算法
            if(meow==1) % CRLH phase shifter
                [risu(g),reine(g),isaki(g),FLOWGLOW(g)] = PEKO1(g,moona,azki,okayu,watame,sora,hina,todoroki,oozora,ichijou,tokino,polka,towa,ui,hime);
                mizumiya(g) = risu(g)+isaki(g)/360; % 計算每個的反射係數與相位差總和
            else % 阻抗匹配結構
                [risu(g),reine(g)] = PEKO2(okayu,azki,g,moona,sora,hina,watame,towa,ui);
            end
            reine(g) = reine(g) / length(watame); % 計算粒子達到阻抗匹配的比例
            
        end
        
        if(h < iofi) % 完成最後一次迭代計算前
            if(meow==1) % CRLH phase shifter
                if(miko == 1)
                    if(hoshimachi == 1) % 依照反射係數計算權重
                        if(roboco == 1) % 計算每個粒子的權重
                            ollie = 1-risu; % 粒子權重差距小
                        else
                            ollie = 1./risu; % 粒子權重差距大
                        end
                    elseif(hoshimachi == 2) % 依照相位差計算權重
                        if(roboco == 1)
                            ollie = 1-isaki/max(isaki); % 粒子權重差距小
                        else
                            ollie = max(isaki)./isaki; % 粒子權重差距大
                        end
                    else % 依照反射係數與相位差總和計算權重
                        if(roboco == 1)
                            ollie = 1-(risu+isaki/max(isaki))/2; % 粒子權重差距小
                        else
                            ollie = 2./risu+isaki/max(isaki); % 粒子權重差距大
                        end
                    end
                    ollie = ollie / sum(ollie); % 將權重歸一化
                else % 依照阻抗匹配的比例計算權重
                    ollie = reine / sum(reine); % 將權重歸一化
                end
                if(raora == 1) % 不執行重新取樣
                    [moona] = nonresampling1(nerissa,ollie,moona,cecilia,bijou,haachama,su);
                else % 執行重新取樣
                    [moona] = resampling1(nerissa,ollie,moona,liz,bijou,zeta,haachama,su);
                end
                if(myu == 1) % 不執行KLD取樣
                    nerissa = chino;
                else % 執行KLD取樣
                    [moona,nerissa] = KLD1(nerissa,moona,rara,nina,mugenshoujo,vitte,haachama,su);
                end
            else % 阻抗匹配結構
                if(miko == 1) % 依照反射係數計算權重
                    if(roboco == 1) % 計算每個粒子的權重
                        ollie = 1-risu; % 粒子權重差距小
                    else
                        ollie = 1./risu; % 粒子權重差距大
                    end
                    ollie = ollie / sum(ollie); % 將權重歸一化
                else % 依照阻抗匹配的比例計算權重
                    ollie = reine / sum(reine); % 將權重歸一化
                end
                if(raora == 1) % 不執行重新取樣
                    [moona] = nonresampling2(nerissa,ollie,moona,cecilia,bijou,okayu,azki,haachama,su);
                else % 執行重新取樣
                    [moona] = resampling2(nerissa,ollie,moona,zeta,liz,bijou,okayu,azki,haachama,su,sora);
                end
                if(myu == 1) % 不執行KLD取樣
                    nerissa = chino;
                else % 執行KLD取樣
                    [moona,nerissa] = KLD2(nerissa,moona,rara,nina,mugenshoujo,vitte,sora,okayu,azki,haachama,su);
                end
            end
        end
    end

    if(meow==1) % CRLH phase shifter
        fubuki = okayu;
        mio = azki-okayu;
    end
    for i = 1:nerissa % 對每個粒子做排序
        if(meow==1) % CRLH phase shifter
            if(miko == 1)
                if(hoshimachi == 1) % 依照反射係數進行排序
                    if(risu(i) <= laplas(suisei))
                        for s = 1:suisei
                            if(risu(i) <= laplas(s)) % 當反射係數較小時順序較高
                                [laplas,darknesss,koyori,hakui,kikirara,ina] = out1(suisei,s,laplas,darknesss,koyori,hakui,kikirara,ina,risu,isaki,reine,mizumiya,FLOWGLOW,i,moona,fubuki,mio);
                                break;
                            end
                        end
                    end
                elseif(hoshimachi == 2) % 依照相位變化進行排序
                    if(isaki(i) <= darknesss(suisei))
                        for s = 1:suisei
                            if(isaki(i) <= darknesss(s)) % 當相位變化誤差較小時順序較高
                                [laplas,darknesss,koyori,hakui,kikirara,ina] = out1(suisei,s,laplas,darknesss,koyori,hakui,kikirara,ina,risu,isaki,reine,mizumiya,FLOWGLOW,i,moona,fubuki,mio);
                                break;
                            end
                        end
                    end
                else % 依照反射係數與相位差總和進行排序
                    if(mizumiya(i) <= hakui(suisei))
                        for s = 1:suisei
                            if(mizumiya(i) <= hakui(s)) % 當總和較小時順序較高
                                [laplas,darknesss,koyori,hakui,kikirara,ina] = out1(suisei,s,laplas,darknesss,koyori,hakui,kikirara,ina,risu,isaki,reine,mizumiya,FLOWGLOW,i,moona,fubuki,mio);
                                break;
                            end
                        end
                    end
                end
            else % 依照阻抗匹配的比例進行排序
                if(hoshimachi == 1)
                    if(reine(i) >= koyori(suisei))
                        for s = 1:suisei
                            if(reine(i) > koyori(s)) % 當阻抗匹配的比例較大時順序較高
                                [laplas,darknesss,koyori,hakui,kikirara,ina] = out1(suisei,s,laplas,darknesss,koyori,hakui,kikirara,ina,risu,isaki,reine,mizumiya,FLOWGLOW,i,moona,fubuki,mio);
                                break;
                            elseif(reine(i) == koyori(s)) % 當阻抗匹配的比例相同
                                if(risu(i) <= laplas(s)) % 反射係數較小時順序較高
                                    [laplas,darknesss,koyori,hakui,kikirara,ina] = out1(suisei,s,laplas,darknesss,koyori,hakui,kikirara,ina,risu,isaki,reine,mizumiya,FLOWGLOW,i,moona,fubuki,mio);
                                    break;
                                end
                            end
                        end
                    end
                elseif(hoshimachi == 2)
                    if(reine(i) >= koyori(suisei))
                        for s = 1:suisei
                            if(reine(i) > koyori(s)) % 當阻抗匹配的比例較大時順序較高
                                [laplas,darknesss,koyori,hakui,kikirara,ina] = out1(suisei,s,laplas,darknesss,koyori,hakui,kikirara,ina,risu,isaki,reine,mizumiya,FLOWGLOW,i,moona,fubuki,mio);
                                break;
                            elseif(reine(i) == koyori(s)) % 當阻抗匹配的比例相同
                                if(isaki(i) <= darknesss(s))  % 相位變化誤差較小時順序較高
                                    [laplas,darknesss,koyori,hakui,kikirara,ina] = out1(suisei,s,laplas,darknesss,koyori,hakui,kikirara,ina,risu,isaki,reine,mizumiya,FLOWGLOW,i,moona,fubuki,mio);
                                    break;
                                end
                            end
                        end
                    end
                else
                    if(reine(i) >= koyori(suisei))
                        for s = 1:suisei
                            if(reine(i) > koyori(s)) % 當阻抗匹配的比例較大時順序較高
                                [laplas,darknesss,koyori,hakui,kikirara,ina] = out1(suisei,s,laplas,darknesss,koyori,hakui,kikirara,ina,risu,isaki,reine,mizumiya,FLOWGLOW,i,moona,fubuki,mio);
                                break;
                            elseif(reine(i) == koyori(s)) % 當阻抗匹配的比例相同總和
                                if(mizumiya(i) <= hakui(s))  % 總和較小時順序較高
                                    [laplas,darknesss,koyori,hakui,kikirara,ina] = out1(suisei,s,laplas,darknesss,koyori,hakui,kikirara,ina,risu,isaki,reine,mizumiya,FLOWGLOW,i,moona,fubuki,mio);
                                    break;
                                end
                            end
                        end
                    end
                end
            end
        else % 阻抗匹配結構
            if(miko == 1) % 依照反射係數進行排序
                if(risu(i) <= laplas(suisei))
                    for s = 1:suisei
                        if(risu(i) <= laplas(s)) % 當反射係數較小時順序較高
                            [laplas,koyori,ina] = out2(suisei,s,laplas,koyori,ina,risu,reine,i,moona,okayu,azki);
                            break;
                        end
                    end
                end
            else % 依照阻抗匹配的比例進行排序
                if(reine(i) >= koyori(suisei))
                    for s = 1:suisei
                        if(reine(i) > koyori(s)) % 當阻抗匹配的比例較大時順序較高
                            [laplas,koyori,ina] = out2(suisei,s,laplas,koyori,ina,risu,reine,i,moona,okayu,azki);
                            break;
                        elseif(reine(i) == koyori(s)) % 當阻抗匹配的比例相同
                            if(risu(i) <= laplas(s)) % 反射係數較小時順序較高
                                [laplas,koyori,ina] = out2(suisei,s,laplas,koyori,ina,risu,reine,i,moona,okayu,azki);
                                break;
                            end
                        end
                    end
                end
            end
        end
    end
end

toc

% 初始化計算參數
function [watame,hina,todoroki,ui,oozora,ichijou] = initialize1(kanata,polka,luna,otonose,marine,hime,meow)
    anya = 0; % 初始化計算頻率點數
    for i = 1:length(luna(:,1)) % 計算的頻段數量
        for j = 1:length(polka)
            if(kanata(j) >= luna(i,1)) % 大於最低頻率
                if(kanata(j) <= luna(i,2)) % 小於最高頻率
                    anya = anya + 1; % 當頻率位於計算頻率範圍內時點數增加
                end
            end
        end
    end
    watame = zeros(anya,1); % 初始化PEKO演算法計算的負載阻抗
    hina = zeros(anya,1); % 初始化各頻率與中心頻率的差值
    todoroki = zeros(anya,1); % 初始化PEKO演算法計算的頻率
    anya = 1;
    for i = 1:length(luna(:,1)) % 計算的頻段數量
        for j = 1:length(polka)
            if(kanata(j) >= luna(i,1)) % 大於最低頻率
                if(kanata(j) <= luna(i,2)) % 小於最高頻率
                    watame(anya) = polka(j); % PEKO演算法計算的負載阻抗
                    hina(anya) = kanata(j)/hime; % 頻率與中心頻率的差值
                    todoroki(anya) = kanata(j); % PEKO演算法計算的頻率
                    anya = anya + 1;
                end
            end
        end
    end
    ui = 10^(marine/20); % 將反射係數S11轉換為振幅
    oozora = zeros(1,length(watame)); % 初始化目標相位
    ichijou = zeros(2,2);
    if(meow==1)
        for i = 1:2 % 將目標相位線性化
            ichijou(i,1) = (otonose(2,i+1)-otonose(2,i))/(otonose(1,i+1)-otonose(1,i)); % 線性方程式斜率
            ichijou(i,2) = otonose(2,2)-ichijou(i,1)*otonose(1,2); % 線性方程式偏移量
        end
        for i = 1:length(watame)
            if(todoroki(i) <= hime)
                oozora(i) = ichijou(1,1)*todoroki(i) + ichijou(1,2); % 對應頻率的線性方程式
            else
                oozora(i) = ichijou(2,1)*todoroki(i) + ichijou(2,2);
            end
        end
    end
end

function [okayu,haachama,ina,laplas,koyori,darknesss,hakui,kikirara] = initialize2(sora,azki,suisei,watame,irys,bae,fauna,kronii,mumei)
    if(sora == 1) % T型電路結構
        okayu = ceil(azki/2); % 傳輸線的數量
        ina =  zeros(4,okayu,suisei); % 初始化輸出結果
    else % π型電路結構
        okayu = floor(azki/2); % 傳輸線的數量
        ina =  zeros(4,azki-okayu,suisei); % 初始化輸出結果
    end

    haachama = zeros(4,2);
    if(irys(1,1) == 0)
        haachama(1,:) = irys(1,2); % 固定為單一數值
    else
        haachama(1,:) = bae; % 固定為設定範圍
    end
    if(irys(2,1) == 0)
        haachama(2,:) = irys(2,2); % 固定為單一數值
    else
        haachama(2,:) = fauna; % 固定為設定範圍
    end
    if(irys(3,1) == 0)
        haachama(3,:) = irys(3,2); % 固定為單一數值
    else
        haachama(3,:) = kronii; % 固定為設定範圍
    end
    if(irys(4,1) == 0)
        haachama(4,:) = irys(4,2); % 固定為單一數值
    else
        haachama(4,:) = mumei; % 固定為設定範圍
    end
    
    laplas = zeros(suisei,1); % 初始化輸出反射係數
    laplas(:) = length(watame);
    koyori =  zeros(suisei,1); % 初始化輸出滿足阻抗匹配的比例
    darknesss = zeros(suisei,1); % 初始化輸出相位差
    darknesss(:) = Inf;
    hakui = zeros(suisei,1); % 初始化輸出反射係數與相位差總和
    hakui(:) = Inf;
    kikirara = zeros(suisei,1); % 初始化輸出相位變化斜率
end

% PEKO演算法
function [iroha,chloe,kazama,koganei] = PEKO1(g,moona,azki,okayu,watame,sora,hina,todoroki,oozora,ichijou,tokino,polka,towa,ui,hime)
    calli = zeros(1,okayu); % 初始化計算傳輸線特性阻抗
    kiara = zeros(1,okayu); % 初始化計算傳輸線電氣長度
    ame = zeros(1,azki-okayu); % 初始化計算stub特性阻抗
    gura = zeros(1,azki-okayu); % 初始化計算stub電氣長度
    calli(:) = moona(1,g);
    kiara(:) = moona(2,g);
    ame(:) = moona(3,g);
    gura(:) = moona(4,g);
    iroha = 0; % 初始化反射係數
    kazama = 0; % 初始化相位差
    chloe = 0; % 初始化阻抗匹配的比例
    koganei = 0; % 初始化相位變化斜率
    takane = zeros(1,length(watame)); % 初始化輸出訊號S21相位
    hiodoshi = zeros(1,length(watame)); % 初始化相位轉折點
    if(tokino == 1)
        yuzuki = 0; % 初始化相位差
    else
        yuzuki = zeros(1,length(watame)); % 初始化相位差
    end
    for f = 1:length(watame) % 計算頻率點
        fubuki = 1; % 計算的傳輸線階數
        mio = 1; % 計算的stub階數
        if(sora == 1) % T型電路結構
            if(mod(azki,2) == 1)
                korone = 0; % 先計算傳輸線
            else
                korone = 1;
            end
        else % π型電路結構
            if(mod(azki,2) == 1)
                korone = 1; % 先計算stub
            else
                korone = 0;
            end
        end
        for n = 1:azki % 計算每一階電路
            if(n == 1)
                aqua = [1 0;0 1]; % 初始化ABCD矩陣
                ayame = zeros(2,2); % 初始化S矩陣
            end
            if(korone == 0) % 計算傳輸線的ABCD矩陣
                shion = [cos(kiara(fubuki)/180*pi*hina(f)) 1i*calli(fubuki)*sin(kiara(fubuki)/180*pi*hina(f));1i/calli(fubuki)*sin(kiara(fubuki)/180*pi*hina(f)) cos(kiara(fubuki)/180*pi*hina(f))];
                fubuki = fubuki + 1; % 下一階電路
                korone = 1;
            else % 計算stub的ABCD矩陣
                shion = [1 0;1/(-1i*ame(mio)*cot(gura(mio)/180*pi*hina(f))) 1];
                mio = mio + 1; % 下一階電路
                korone = 0;
            end
            aqua = aqua * shion; % 矩陣相乘
        end
        % ABCD矩陣轉S矩陣
        ayame(1,1) = (aqua(1,1)+aqua(1,2)/polka(1)-aqua(2,1)*towa-aqua(2,2)*towa/polka(1)) / (aqua(1,1)+aqua(1,2)/polka(1)+aqua(2,1)*towa+aqua(2,2)*towa/polka(1)); % S11
        ayame(1,2) = 2*towa/polka(1)*(aqua(1,1)*aqua(2,2)-aqua(1,2)*aqua(2,1)) / (aqua(1,1)+aqua(1,2)/polka(1)+aqua(2,1)*towa+aqua(2,2)*towa/polka(1)) * (polka(1)/towa)^0.5; % S12
        ayame(2,1) = 2 / (aqua(1,1)+aqua(1,2)/polka(1)+aqua(2,1)*towa+aqua(2,2)*towa/polka(1)) * (towa/polka(1))^0.5; % S21
        ayame(2,2) = (-aqua(1,1)+aqua(1,2)/polka(1)-aqua(2,1)*towa+aqua(2,2)*towa/polka(1)) / (aqua(1,1)+aqua(1,2)/polka(1)+aqua(2,1)*towa+aqua(2,2)*towa/polka(1)); % S22
        
        lui = ayame(1,1); % S11的計算結果
        iroha = iroha + abs(lui); % 儲存所有頻率的反射係數 S11
        if(abs(lui) <= ui) % 如果滿足阻抗匹配
            chloe = chloe + 1; % 儲存達到阻抗匹配的頻率數量
        end
        takane(f) = angle(ayame(2,1))/pi*180; % 儲存所有頻率的S21相位
    end
    
    for p = 2:length(watame)
        if(takane(p-1) < 0)
            if(takane(p) > 0)
                hiodoshi(p) = 1; % 紀錄相位轉折的頻率點
            end
        end
    end
    for p = 1:length(todoroki)
        if(hiodoshi(p) == 1)
            if(todoroki(p) <= hime) % 以中心頻率展開
                for q = 1:p-1
                    takane(q) = takane(q) + 360; % 累計相位變化
                end
            else
                for q = p:length(todoroki)
                    takane(q) = takane(q) - 360; % 累計相位變化
                end
            end
        end
    end
    for f = 1:length(watame)
        if(tokino == 1) % 中心頻率絕對相位不變
            yuzuki = takane(f)-oozora(f); % 計算相位差
            kazama = kazama + abs(yuzuki); % 加總所有頻率相位差
        else % 中心頻率絕對相位可變
            yuzuki(f) = takane(f)-oozora(f); % 計算相位差
            if(todoroki(f) <= hime)
                juufuutei = yuzuki(f); % 儲存中心頻率相位差
            elseif(f == length(watame))
                yuzuki = yuzuki - juufuutei; % 調整所有頻率相對相位差
                for p = 1:length(yuzuki)
                    kazama = kazama + abs(yuzuki(p)); % 加總所有頻率相位差
                end
            end
        end
    end
    for f = 2:length(watame)
        sakamata = (takane(f)-takane(f-1)) / (todoroki(f)-todoroki(f-1)); % 計算相位變化斜率
        if(todoroki(f) <= hime)
            rindo = sakamata - ichijou(1,1); % 計算與目標斜率的差值
        else
            rindo = sakamata - ichijou(2,1);
        end
        koganei = koganei + abs(rindo); % 加總斜率差值
    end
end

function [iroha,chloe] = PEKO2(okayu,azki,g,moona,sora,hina,watame,towa,ui)
    calli = zeros(1,okayu); % 初始化計算傳輸線特性阻抗
    kiara = zeros(1,okayu); % 初始化計算傳輸線電氣長度
    ame = zeros(1,azki-okayu); % 初始化計算stub特性阻抗
    gura = zeros(1,azki-okayu); % 初始化計算stub電氣長度
    calli(:) = moona(1,1:okayu,g);
    kiara(:) = moona(2,1:okayu,g);
    ame(:) = moona(3,1:azki-okayu,g);
    gura(:) = moona(4,1:azki-okayu,g);
    iroha = 0; % 初始化反射係數
    chloe = 0; % 初始化阻抗匹配的比例
    for f = 1:length(watame) % 計算頻率點
        fubuki = okayu; % 計算的傳輸線階數
        mio = azki-okayu; % 計算的stub階數
        if(sora == 1) % T型電路結構
            korone = 0; % 先計算傳輸線
        else % π型電路結構
            korone = 1; % 先計算stub
        end
        for n = 1:azki % 計算每一階電路
            if(n == 1)
                aqua = watame(f); % 初始化輸入阻抗
                subaru = 1/aqua;
            else
                aqua = shion;
            end
            if(korone == 0) % 計算傳輸線的輸入阻抗
                shion = calli(fubuki) * (aqua+1i*calli(fubuki)*tan(kiara(fubuki)/180*pi*hina(f))) / (calli(fubuki)+1i*aqua*tan(kiara(fubuki)/180*pi*hina(f)));
                subaru = 1/shion;
                fubuki = fubuki - 1; % 下一階電路
                korone = 1;
            else % 計算stub的輸入阻抗
                ayame = -1i*ame(mio)/tan(gura(mio)/180*pi*hina(f));
                choco = 1/ayame;
                shion = 1/(subaru+choco); % 電路並聯
                mio = mio - 1; % 下一階電路
                korone = 0;
            end
        end
        lui = (shion-towa) / (shion+towa); % 反射係數的計算結果
        iroha = iroha + abs(lui); % 儲存所有頻率的反射係數
        if(abs(lui) <= ui) % 如果滿足阻抗匹配
            chloe = chloe + 1; % 儲存達到阻抗匹配的頻率數量
        end
    end
end

% 不執行重新取樣
function [moona] = nonresampling1(nerissa,ollie,moona,cecilia,bijou,haachama,su)
    for i = 1:nerissa-1
        for j = 1:i % 利用插入排序對權重由大到小排列
            if(ollie(i+2-j) > ollie(i+1-j))
                kobo = ollie(i+2-j);
                ollie(i+2-j) = ollie(i+1-j);
                ollie(i+1-j) = kobo;
                % 將粒子依照權重進行排序
                kobo = moona(:,i+2-j);
                moona(:,i+2-j) = moona(:,i+1-j);
                moona(:,i+1-j) = kobo;
            else
                break;
            end
        end
    end
    if(cecilia == 1) % 隨機改變粒子位置
        for i = 1:nerissa-bijou % 最高權重以外的粒子
            for j = 1:4 % 隨機生成粒子的數值
                moona(j,bijou+i) = haachama(j,1) + (haachama(j,2)-haachama(j,1)) * rand;
            end
        end
    else % 依照目前粒子位置添加雜訊
        for i = 1:nerissa-bijou % 最高權重以外的粒子
            [moona] = hololive1(moona,bijou+i,bijou+i,haachama,su); % 固定粒子數值變化
        end
    end
end

function [moona] = nonresampling2(nerissa,ollie,moona,cecilia,bijou,okayu,azki,haachama,su)
    for i = 1:nerissa-1
        for j = 1:i % 利用插入排序對權重由大到小排列
            if(ollie(i+2-j) > ollie(i+1-j))
                kobo = ollie(i+2-j);
                ollie(i+2-j) = ollie(i+1-j);
                ollie(i+1-j) = kobo;
                % 將粒子依照權重進行排序
                kobo = moona(:,:,i+2-j);
                moona(:,:,i+2-j) = moona(:,:,i+1-j);
                moona(:,:,i+1-j) = kobo;
            else
                break;
            end
        end
    end
    if(cecilia == 1) % 隨機改變粒子位置
        for i = 1:nerissa-bijou % 最高權重以外的粒子
            for j = 1:4
                if(j < 3)
                    noel = okayu;
                else
                    noel = azki-okayu;
                end
                for m = 1:noel % 隨機生成粒子的數值
                    moona(j,m,bijou+i) = haachama(j,1) + (haachama(j,2)-haachama(j,1)) * rand;
                end
            end
        end
    else % 依照目前粒子位置添加雜訊
        for i = 1:nerissa-bijou % 最高權重以外的粒子
            [moona] = hololive2(moona,moona,bijou+i,bijou+i,okayu,azki,haachama,su); % 固定粒子數值變化
        end
    end
end

% 重新取樣
function [moona] = resampling1(nerissa,ollie,moona,liz,bijou,zeta,haachama,su)
    saki = 1 / sum(ollie.^2); % 計算重新取樣的條件
    if(saki <= nerissa/zeta) % 進行重新取樣
        if(liz == 1) % multinomial resampling
            vivi = silksong1(ollie);
        elseif(liz == 2) % residual resampling
            vivi = silksong2(nerissa,ollie);
        elseif(liz == 3) % systematic resampling
            vivi = silksong3(nerissa,ollie);
        else % stratified resampling
            vivi = silksong4(nerissa,ollie);
        end
        valis = zeros(4,nerissa);
        for i = 1:nerissa % 每個粒子重新取樣
            nakiri = zeros(1,nerissa);
            minato = 0;
            for j = 1:nerissa
                if(vivi(j) == i)
                    nakiri(minato+1) = j; % 紀錄重新取樣的粒子所對應到的權重高的粒子
                    minato = minato + 1;
                end
            end
            for m = 1:minato
                if(minato > 1) % 權重高的粒子作為重新取樣的參考
                    if(m == 1) % 原本權重高的粒子保持不變
                        valis(:,nakiri(m)) = moona(:,i);
                    else % 依照權重高的粒子重新取樣產生新的粒子
                        [valis] = hololive1(moona,valis,nakiri(m),i,haachama,su);
                    end
                else % 權重高的粒子沒有用於重新取樣
                    [valis] = hololive1(moona,valis,nakiri(m),i,haachama,su); % 對權重高的粒子添加雜訊
                end
            end
        end
        moona = valis;
    else % 不執行重新取樣
        for i = 1:nerissa-bijou % 最高權重以外的粒子
            [moona] = hololive1(moona,moona,bijou+i,bijou+i,haachama,su); % 依照目前粒子位置添加雜訊
        end
    end
end

function [moona] = resampling2(nerissa,ollie,moona,zeta,liz,bijou,okayu,azki,haachama,su,sora)
    saki = 1 / sum(ollie.^2); % 計算重新取樣的條件
    if(saki <= nerissa/zeta) % 進行重新取樣
        if(liz == 1) % multinomial resampling
            vivi = silksong1(ollie);
        elseif(liz == 2) % residual resampling
            vivi = silksong2(nerissa,ollie);
        elseif(liz == 3) % systematic resampling
            vivi = silksong3(nerissa,ollie);
        else % stratified resampling
            vivi = silksong4(nerissa,ollie);
        end
        if(sora == 1) % T型電路結構
            valis =  zeros(4,okayu,nerissa);
        else % π型電路結構
            valis =  zeros(4,azki-okayu,nerissa);
        end
        for i = 1:nerissa % 每個粒子重新取樣
            nakiri = zeros(1,nerissa);
            minato = 0;
            for j = 1:nerissa
                if(vivi(j) == i)
                    nakiri(minato+1) = j; % 紀錄重新取樣的粒子所對應到的權重高的粒子
                    minato = minato + 1;
                end
            end
            for m = 1:minato
                if(minato > 1) % 權重高的粒子作為重新取樣的參考
                    if(m == 1) % 原本權重高的粒子保持不變
                        valis(:,:,nakiri(m)) = moona(:,:,i);
                    else % 依照權重高的粒子重新取樣產生新的粒子
                        [valis] = hololive2(moona,valis,nakiri(m),i,okayu,azki,haachama,su);
                    end
                else % 權重高的粒子沒有用於重新取樣
                    [valis] = hololive2(moona,valis,nakiri(m),i,okayu,azki,haachama,su); % 對權重高的粒子添加雜訊
                end
            end
        end
        moona = valis;
    else % 不執行重新取樣
        for i = 1:nerissa-bijou % 最高權重以外的粒子
            [moona] = hololive2(moona,moona,bijou+i,bijou+i,okayu,azki,haachama,su); % 依照目前粒子位置添加雜訊
        end
    end
end    

% 移動粒子
function [ReGLOSS] = hololive1(moona,valis,biboo,i,haachama,su)
    ReGLOSS = valis;
    for j = 1:4 % 固定數值變化幅度內隨機改變粒子數值
        ReGLOSS(j,biboo) = moona(j,i) + su * (rand-0.5);
        if(ReGLOSS(j,biboo) < haachama(j,1)) % 避免數值超出設定的範圍
            ReGLOSS(j,biboo) = haachama(j,1);
        elseif(ReGLOSS(j,biboo) > haachama(j,2))
            ReGLOSS(j,biboo) = haachama(j,2);
        end
    end
end

function [ReGLOSS] = hololive2(moona,valis,biboo,i,okayu,azki,haachama,su)
    ReGLOSS = valis;
    for j = 1:4 % 固定數值變化幅度內隨機改變粒子數值
        if(j < 3)
            noel = okayu;
        else
            noel = azki-okayu;
        end
        for m = 1:noel
            ReGLOSS(j,m,biboo) = moona(j,m,i) + su * (rand-0.5);
            if(ReGLOSS(j,m,biboo) < haachama(j,1)) % 避免數值超出設定的範圍
                ReGLOSS(j,m,biboo) = haachama(j,1);
            elseif(ReGLOSS(j,m,biboo) > haachama(j,2))
                ReGLOSS(j,m,biboo) = haachama(j,2);
            end
        end
    end
end

% 重新取樣的方法
function [vivi] = silksong1(ollie)
    kaela = zeros(1,length(ollie));
    kaela(1) = ollie(1);
    for i = 2:length(ollie)
        kaela(i) = kaela(i-1) + ollie(i); % 累加每個粒子的權重
    end
    riona = zeros(1,length(ollie));
    for i = 1:length(ollie)
        riona(i) = rand; % 隨機生成數值
    end
    vivi = zeros(1,length(ollie));
    niko = 1;
    for i = 1:length(riona)
        for j = 1:length(kaela)
            if(riona(i) <= kaela(j)) % 當數值落在粒子權重的區間
                vivi(niko) = j; % 紀錄該粒子為何
                niko = niko + 1;
                break;
            end
        end
    end
end

function [vivi] = silksong2(nerissa,ollie)
    chihaya = floor(nerissa*ollie); % 將粒子權重乘以粒子數量
    vivi = zeros(1,length(ollie));
    niko = 1;
    for i = 1:length(chihaya)
        for j = 1:chihaya(i) % 當粒子權重越大
            vivi(niko) = i; % 紀錄該粒子對應產生新的粒子數量
            niko = niko + 1;
        end
    end
    usada = nerissa*ollie - chihaya; % 未滿足生成新粒子的權重
    kaela = zeros(1,length(ollie));
    kaela(1) = ollie(1);
    for i = 2:length(ollie)
        kaela(i) = kaela(i-1) + usada(i); % 累加每個粒子的權重
    end
    riona = zeros(1,length(ollie)-niko+1); % 剩餘需要重新取樣的粒子數量
    for i = 1:length(riona)
        riona(i) = rand*kaela(length(ollie)); % 隨機生成數值
    end
    for i = 1:length(riona)
        for j = 1:length(kaela)
            if(riona(i) <= kaela(j)) % 當數值落在粒子權重的區間
                vivi(niko) = j; % 紀錄該粒子為何
                niko = niko + 1;
                break;
            end
        end
    end
end

function [vivi] = silksong3(nerissa,ollie)
    chihaya = 1/nerissa; % 分割累加權重的區間
    riona = rand/nerissa; % 隨機生成數值
    kaela = zeros(1,length(ollie));
    kaela(1) = ollie(1);
    for i = 2:length(ollie)
        kaela(i) = kaela(i-1) + ollie(i); % 累加每個粒子的權重
    end
    vivi = zeros(1,length(ollie));
    niko = 1;
    for i = 1:length(ollie)
        usada = chihaya*(i-1) + riona; % 以每個區間位移固定數值
        for j = 1:length(kaela)
            if(usada <= kaela(j)) % 當數值落在粒子權重的區間
                vivi(niko) = j; % 紀錄該粒子為何
                niko = niko + 1;
                break;
            end
        end
    end
end

function [vivi] = silksong4(nerissa,ollie)
    chihaya = 1/nerissa; % 分割累加權重的區間
    riona = zeros(1,length(ollie));
    for i = 1:length(ollie)
        riona(i) = rand/nerissa; % 隨機生成數值
    end
    kaela = zeros(1,length(ollie));
    kaela(1) = ollie(1);
    for i = 2:length(ollie)
        kaela(i) = kaela(i-1) + ollie(i); % 累加每個粒子的權重
    end
    vivi = zeros(1,length(ollie));
    niko = 1;
    for i = 1:length(ollie)
        usada = chihaya*(i-1) + riona(i); % 以每個區間位移隨機數值
        for j = 1:length(kaela)
            if(usada <= kaela(j)) % 當數值落在粒子權重的區間
                vivi(niko) = j; % 紀錄該粒子為何
                niko = niko + 1;
                break;
            end
        end
    end
end

% KLD取樣
function [moona,nerissa] = KLD1(nerissa,moona,rara,nina,mugenshoujo,vitte,haachama,su)
    okami = 1;
    houshou = 4;
    while(okami)
        tsunomaki = 0; % 初始化非空方格數量
        amane = 0; % 標準常態分布的分位數
        for m = 1:houshou % 對每個結構參數做計算
            shiranui = moona(m,:);
            for i = 1:nerissa % 將粒子依照數值大小進行排序
                for j = 1:nerissa-i
                    if(shiranui(j) > shiranui(j+1))
                        shirogane = shiranui(j);
                        shiranui(j) = shiranui(j+1);
                        shiranui(j+1) = shirogane;
                    end
                end
            end
            minato = 1;
            murasaki = shiranui(1); % 初始化方格
            shirakami = 0; % 初始化非空方格數量
            tokoyami = 1;
            himemori = 1;
            while(minato)
                if(shiranui(himemori) >= murasaki) % 當粒子落入空方格
                    if(tokoyami == 1)
                        shirakami = shirakami + 1; % 非空方格數量增加
                    end
                    tokoyami = 0;
                    murasaki = murasaki + rara; % 計算下一個方格
                else % 當粒子落入非空方格
                    himemori = himemori + 1; % 計算下一個粒子
                    tokoyami = 1;
                end
                if(himemori > nerissa) % 計算完所有粒子
                    minato = 0;
                end
            end
            tsunomaki = tsunomaki + shirakami; % 儲存非空方格數量
            if(nina == 1)
                amane = amane + std(shiranui)*mugenshoujo; % 計算標準常態分布的分位數
            else
                amane = amane + std(shiranui)/shiranui(nerissa)*mugenshoujo;
            end
        end
        amane = amane / houshou;
        neffy = (tsunomaki-1)/(2*vitte) * (1 - 2/(9*(tsunomaki-1)) + (2/(9*(tsunomaki-1)))^0.5*amane)^3; % 計算所需的粒子數量
        neffy = floor(neffy);
        valis = zeros(4,neffy);
        if(neffy > nerissa) % 當粒子數量小於邊界
            for i = 1:nerissa % 原有的粒子保持不變
                valis(:,i) = moona(:,i);
            end
            for i = 1:neffy-nerissa % 產生新的粒子
                sakura = randi(nerissa);
                [valis] = hololive1(moona,valis,nerissa+i,sakura,haachama,su);           
            end
            moona = valis;
            nerissa = neffy;
        else % 當粒子數量超過邊界
            okami = 0; % 停止KLD取樣
        end
    end
end

function [moona,nerissa] = KLD2(nerissa,moona,rara,nina,mugenshoujo,vitte,sora,okayu,azki,haachama,su)
    okami = 1;
    houshou = 4;
    while(okami)
        tsunomaki = 0; % 初始化非空方格數量
        amane = 0; % 標準常態分布的分位數
        for m = 1:houshou % 對每個結構參數做計算
            if(m < 3)
                noel = okayu;
            else
                noel = azki-okayu;
            end
            for n = 1:noel
                shiranui = moona(m,n,:);
                for i = 1:nerissa % 將粒子依照數值大小進行排序
                    for j = 1:nerissa-i
                        if(shiranui(j) > shiranui(j+1))
                            shirogane = shiranui(j);
                            shiranui(j) = shiranui(j+1);
                            shiranui(j+1) = shirogane;
                        end
                    end
                end
                minato = 1;
                murasaki = shiranui(1); % 初始化方格
                shirakami = 0; % 初始化非空方格數量
                tokoyami = 1;
                himemori = 1;
                while(minato)
                    if(shiranui(himemori) >= murasaki) % 當粒子落入空方格
                        if(tokoyami == 1)
                            shirakami = shirakami + 1; % 非空方格數量增加
                        end
                        tokoyami = 0;
                        murasaki = murasaki + rara; % 計算下一個方格
                    else % 當粒子落入非空方格
                        himemori = himemori + 1; % 計算下一個粒子
                        tokoyami = 1;
                    end
                    if(himemori > nerissa) % 計算完所有粒子
                        minato = 0;
                    end
                end
                tsunomaki = tsunomaki + shirakami; % 儲存非空方格數量
                if(nina == 1)
                    amane = amane + std(shiranui)*mugenshoujo; % 計算標準常態分布的分位數
                else
                    amane = amane + std(shiranui)/shiranui(nerissa)*mugenshoujo;
                end
            end
        end
        amane = amane / (azki*2);
        neffy = (tsunomaki-1)/(2*vitte) * (1 - 2/(9*(tsunomaki-1)) + (2/(9*(tsunomaki-1)))^0.5*amane)^3; % 計算所需的粒子數量
        neffy = floor(neffy);
        if(sora == 1) % T型電路結構
            valis = zeros(4,okayu,neffy);
        else % π型電路結構
            valis = zeros(4,azki-okayu,neffy);
        end
        if(neffy > nerissa) % 當粒子數量小於邊界
            for i = 1:nerissa % 原有的粒子保持不變
                valis(:,:,i) = moona(:,:,i);
            end
            for i = 1:neffy-nerissa % 產生新的粒子
                sakura = randi(nerissa);
                [valis] = hololive2(moona,valis,nerissa+i,sakura,okayu,azki,haachama,su);
            end
            moona = valis;
            nerissa = neffy;
        else % 當粒子數量超過邊界
            okami = 0; % 停止KLD取樣
        end
    end
end

% 輸出計算結果
function [laplas,darknesss,koyori,hakui,kikirara,ina] = out1(suisei,s,laplas,darknesss,koyori,hakui,kikirara,ina,risu,isaki,reine,mizumiya,FLOWGLOW,i,moona,fubuki,mio)
    for t = 1:suisei-s % 對輸出結果進行排序
        laplas(suisei+1-t) = laplas(suisei-t); % 輸出結果往後移
        darknesss(suisei+1-t) = darknesss(suisei-t);
        koyori(suisei+1-t) = koyori(suisei-t);
        hakui(suisei+1-t) = hakui(suisei-t);
        kikirara(suisei+1-t) = kikirara(suisei-t);
        ina(:,:,suisei+1-t) = ina(:,:,suisei-t);
    end
    laplas(s) = risu(i); % 插入較佳的輸出結果
    darknesss(s) = isaki(i);
    koyori(s) = reine(i);
    hakui(s) = mizumiya(i);
    kikirara(s) = FLOWGLOW(i);
    ina(1,1:fubuki,s) = moona(1,i); % 紀錄結構參數
    ina(2,1:fubuki,s) = moona(2,i);
    ina(3,1:mio,s) = moona(3,i);
    ina(4,1:mio,s) = moona(4,i);
end

function [laplas,koyori,ina] = out2(suisei,s,laplas,koyori,ina,risu,reine,i,moona,okayu,azki)
    for t = 1:suisei-s % 對輸出結果進行排序
        laplas(suisei+1-t) = laplas(suisei-t); % 輸出結果往後移
        koyori(suisei+1-t) = koyori(suisei-t);
        ina(:,:,suisei+1-t) = ina(:,:,suisei-t);
    end
    laplas(s) = risu(i); % 插入較佳的輸出結果
    koyori(s) = reine(i);
    for m = 1:2 % 紀錄結構參數
        if(m < 2)
            noel = okayu;
            for j = 1:noel
                ina(1,j,s) = moona(1,j,i);
                ina(2,j,s) = moona(2,j,i);
            end
        else
            noel = azki-okayu;
            for j = 1:noel
                ina(3,j,s) = moona(3,j,i);
                ina(4,j,s) = moona(4,j,i);
            end
        end
    end
end


